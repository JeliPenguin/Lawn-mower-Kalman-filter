function ins_gnss_solution = INS_GNSS(gnss_solution,dr_gnss_solution,dr_measurement_data,times)
% Calculates the integrated solution of GNSS and Inertial Navigation System
% Inputs:
%   gnss_solution                             Calculated GNSS Kalman Filter solution
%   dr_solution                               GNSS corrected Dead Reckoning solution
%   times                                     Array of time value for each epoch
%
% Outputs:
%   ins_gnss_solution                         Calculated DR/GNSS solution
%                                               Col 1: Time (s)
%                                               Col 2: Latitude (deg)
%                                               Col 3: Longitude (deg)
%                                               Col 4: North Velocity (m/s)
%                                               Col 5: East Velocity (m/s)
%                                               Col 6: Heading (deg)

Define_Constants;

dr_measurement = dr_measurement_data(:,2:end);
dr_gnss_solution = table2array(inertiaNavigationTable);
epochNum = height(times);

% Initialise state estimates
x_caret_k_m_1_plus = zeros(15,1);

% Initialise error covariance matrix
P_k_m_1_plus =  zeros(15);
P_k_m_1_plus(1:3,1:3) = eye(3) * deg_to_rad^2;
P_k_m_1_plus(4:6,4:6) = eye(3) * 0.1^2;
P_k_m_1_plus(7:9,7:9) = eye(3) * 10^2;
P_k_m_1_plus(10:12,10:12) = eye(3) *...
    (1000 * micro_g_to_meters_per_second_squared)^2;
P_k_m_1_plus(13:15,13:15) = eye(3) * (10 * deg_to_rad / 3600)^2;

ins_gnss_solution = [];

Q_k_m_1 = [LC_KF_config.gyro_noise_PSD*tau_s*eye(3),zeros(3),zeros(3),zeros(3),zeros(3)
           zeros(3),LC_KF_config.accel_noise_PSD*tau_s*eye(3),zeros(3),zeros(3),zeros(3)
           zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)
           zeros(3),zeros(3),zeros(3),LC_KF_config.accel_bias_PSD*tau_s*eye(3),zeros(3)
           zeros(3),zeros(3),zeros(3),zeros(3),LC_KF_config.gyro_bias_PSD*tau_s*eye(3)
           ];


for i=1:epochNum

    time = times(i);
    L_b = dr_gnss_solution(i,2) * deg_to_rad;
    lambda_b = dr_gnss_solution(i,3) * deg_to_rad;
    h_b = gnss_solution(i,4); % No height readings from DR file

    v_eb_n = [dr_gnss_solution(i,5:6)';0]; % Adding arbitrary v_D as we are computing horizontal solutions

    f_ib_b = dr_gnss_solution(i,11:end)'; 

    gnss_lat = gnss_solution(i,2)*deg_to_rad;
    gnss_long = gnss_solution(i,3)*deg_to_rad;
    gnss_h = gnss_solution(i,4);
    gnss_v = gnss_solution(i,5:7)';

    % roll,pitch,heading
    heading = dr_measurement(i,6) * deg_to_rad;
    gyro = [0,0,heading];  % Land vehicle thus assuming approx. 0 roll and pitch

    C_n_b = Euler_to_CTM(gyro);
    C_b_n = C_n_b';

    [r_eb_e,v_eb_e,C_b_e] = NED_to_ECEF(L_b,lambda_b,h_b,v_eb_n,C_b_n);

    F_21 = Calculate_F21(C_b_e,f_ib_b);
    F_23 = Calculate_F23(r_eb_e,L_b);

    Phi_k_m_1 = [eye(3)-Omega_ie*tau_s,zeros(3),zeros(3),zeros(3),C_b_e*tau_s
                 F_21*tau_s,eye(3)-2*Omega_ie*tau_s,F_23*tau_s,C_b_e*tau_s,zeros(3)
                 zeros(3),eye(3)*tau_s,eye(3),zeros(3),zeros(3)
                 zeros(3),zeros(3),zeros(3),eye(3),zeros(3)
                 zeros(3),zeros(3),zeros(3),zeros(3),eye(3)
                 ];

    x_caret_k_minus = Phi_k_m_1*x_caret_k_m_1_plus;
    P_k_minus = Phi_k_m_1*P_k_m_1_plus*Phi_k_m_1'+Q_k_m_1;

    H_k = [zeros(3),zeros(3),-eye(3),zeros(3),zeros(3)
           zeros(3),-eye(3),zeros(3),zeros(3),zeros(3)];

    R_k = [LC_KF_config.pos_meas_SD^2*eye(3),zeros(3)
           zeros(3),LC_KF_config.vel_meas_SD^2*eye(3)
           ];

    K_k = P_k_minus*H_k'*inv(H_k*P_k_minus*H_k'+R_k);

    [r_caret_ea_e,v_caret_ea_e] = pv_NED_to_ECEF(gnss_lat,gnss_long,gnss_h,gnss_v);

    delta_z_k_minus = [r_caret_ea_e;v_caret_ea_e] - [r_eb_e;v_eb_e] + [x_caret_k_minus(7:9);x_caret_k_minus(4:6)];

    x_caret_k_m_1_plus = x_caret_k_minus + K_k*delta_z_k_minus;
    P_k_m_1_plus = (eye(15)-K_k*H_k)*P_k_minus;

    skewed = Skew_symmetric(x_caret_k_m_1_plus(1:3));
    C_caret_b_e = (eye(3)-skewed)*C_b_e;
    v_caret_eb_e = v_eb_e - x_caret_k_m_1_plus(4:6);
    r_caret_eb_e = r_eb_e - x_caret_k_m_1_plus(7:9);

    [L_b,lambda_b,h_b,v_eb_n,C_b_n] = ECEF_to_NED(r_caret_eb_e,v_caret_eb_e,C_caret_b_e);
    C_n_b = C_b_n.';
    eul = CTM_to_Euler(C_n_b)*rad_to_deg;

    ins_gnss_solution = [ins_gnss_solution;[time,L_b*rad_to_deg,lambda_b*rad_to_deg,v_eb_n(1:2)',eul']];

end


% Write to table and export
outputTable = table(array2table(ins_gnss_solution));
writetable(outputTable,"Solutions/INS_GNSS_Solution.csv",'WriteVariableNames',0)
end
