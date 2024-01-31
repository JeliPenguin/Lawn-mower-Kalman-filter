function gyro_mag_heading_solution = Gyro_Magnetometer_KF(dr_measurement_data,times)

Define_Constants;

dr_measurement = dr_measurement_data(:,2:end);

epoch_num = height(times);

P_k_m_1_plus = [sigma_v^2,0,0,0
                0,sigma_v^2,0,0
                0,0,sigma_r^2/(R_N+gnss_h)^2,0
                0,0,0,sigma_r^2/((R_E+gnss_h)^2*(cos(gnss_lat))^2)];


x_caret_k_m_1_plus = zeros(2,1);

% Compute transition matrix
Phi_k_m_1 = [1,tau
             0,1];

% Compute system noise covariance matrix
Q_k_m_1 = [S_rg*tau+(S_bgd*tau^3)/3,0.5*S_bgd*tau^2
           0.5*S_bgd*tau^2,S_bgd*tau
           ];

gyro_mag_heading_solution = zeros(epoch_num,1);
gyro_mag_heading_solution(1) = dr_measurement(1:7);

for i = 2:epoch_num

    time = times(i);

    magnetic_heading_measurement = dr_measurement(i:7);
    gyro_measurement = dr_measurement(i:6);

    % Propagate state estimates and error covariance matrix
    x_caret_k_minus = Phi_k_m_1*x_caret_k_m_1_plus;
    P_k_minus = Phi_k_m_1*P_k_m_1_plus*Phi_k_m_1'+Q_k_m_1;

    % Compute measurement matrix 
    H_k = [-1,0];

    % Compute measurement noise covariance
    R_k = sigma_mh^2;

    % Compute Kalman gain matrix
    K_k = P_k_minus*H_k'*inv(H_k*P_k_minus*H_k'+R_k);

    % Formulate measurement innovation vector
    delta_z_m_k = psi_k_M-psi_k_g+x_caret_k_minus(1);

    % Update state estimates and error covaraince matrix
    x_caret_k_m_1_plus = x_caret_k_minus + K_k*delta_z_m_k;
    P_k_m_1_plus = (eye(2) - K_k*H_k)*P_k_minus;

    % Apply correction to solution of current epoch
    correction = []-x_caret_k_m_1_plus;
    gyro_mag_heading_solution(i) = correction;

end

% Writing solution to csv file
outputTable = array2table(dr_gnss_solution);
writetable(outputTable,"Solutions/DR_GNSS_Solution.csv",'WriteVariableNames',0)

end