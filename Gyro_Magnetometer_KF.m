function gyro_mag_heading_solution = Gyro_Magnetometer_KF(dr_measurement_data,times)

Define_Constants;

dr_measurement = dr_measurement_data(:,2:end);

epoch_num = height(times);

P_k_m_1_plus = [sigma_gyro^2,0
                0,1];


x_caret_k_m_1_plus = zeros(2,1);

% Compute transition matrix
Phi_k_m_1 = [1,tau
             0,1];

% Compute system noise covariance matrix
Q_k_m_1 = [S_rg*tau+(S_bgd*tau^3)/3,0.5*S_bgd*tau^2
           0.5*S_bgd*tau^2,S_bgd*tau
           ];

% Compute measurement matrix 
H_k = [-1,0];

% Compute measurement noise covariance
R_k = sigma_mh^2;

gyro_data = dr_measurement(:,5); 
mag_data = dr_measurement(:,6)*deg_to_rad;

gyro_mag_heading_solution = zeros(epoch_num,1);
gyro_mag_heading_solution(1) = mag_data(1)*rad_to_deg;
psi_k_g = mag_data(1);
for i = 2:epoch_num

    psi_k_M = mag_data(i);
    psi_k_g = psi_k_g+gyro_data(i)*tau;

    % Propagate state estimates and error covariance matrix
    x_caret_k_minus = Phi_k_m_1*x_caret_k_m_1_plus;
    P_k_minus = Phi_k_m_1*P_k_m_1_plus*Phi_k_m_1'+Q_k_m_1;

    % Compute Kalman gain matrix
    K_k = P_k_minus*H_k'*inv(H_k*P_k_minus*H_k'+R_k);

    % Formulate measurement innovation vector
    delta_z_m_k = psi_k_M-psi_k_g+x_caret_k_minus(1);

    % Update state estimates and error covaraince matrix
    x_caret_k_m_1_plus = x_caret_k_minus + K_k*delta_z_m_k;
    P_k_m_1_plus = (eye(2) - K_k*H_k)*P_k_minus;

    % Apply correction to solution of current epoch
    correction = psi_k_g-x_caret_k_m_1_plus(1);
    gyro_mag_heading_solution(i) = correction*rad_to_deg;

end

t = (1:epoch_num) * tau;
f = figure("Visible","off");
plot(t, mag_data*rad_to_deg, 'b', t, gyro_mag_heading_solution, 'r');
legend('Magnetometer Data', 'Filtered Heading');
xlabel('Time (s)');
ylabel('Heading (deg)');
title('Gyro Smoothed magnetic heading');
saveas(f,"Figures/Heading/smoothedComparison.png")

% Writing solution to csv file
outputTable = array2table(gyro_mag_heading_solution);
writetable(outputTable,"Solutions/gyro_smoothed_heading_solution.csv",'WriteVariableNames',0)

end