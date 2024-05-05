clear;
Define_Constants;

task2SolutionTable = readtable("WS2_Task2Solution.csv");
gnssPosVel = table2array(task2SolutionTable);

inertiaNavigationTable = readtable("Aircraft_inertial.csv");
inertialNavigation = table2array(inertiaNavigationTable);
epochNum = height(inertialNavigation);

[x_caret_k_m_1_plus,P_k_m_1_plus] = Initialise_Integration_KF;


% Define Constants
tau_s = 1;

solution = [];

Q_k_m_1 = [LC_KF_config.gyro_noise_PSD*tau_s*eye(3),zeros(3),zeros(3),zeros(3),zeros(3)
           zeros(3),LC_KF_config.accel_noise_PSD*tau_s*eye(3),zeros(3),zeros(3),zeros(3)
           zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)
           zeros(3),zeros(3),zeros(3),LC_KF_config.accel_bias_PSD*tau_s*eye(3),zeros(3)
           zeros(3),zeros(3),zeros(3),zeros(3),LC_KF_config.gyro_bias_PSD*tau_s*eye(3)
           ];


for i=1:epochNum

    time = inertialNavigation(i,1);
    L_b = inertialNavigation(i,2) * deg_to_rad;
    lambda_b = inertialNavigation(i,3) * deg_to_rad;
    h_b = inertialNavigation(i,4);
    v_eb_n = inertialNavigation(i,5:7)';
    gyro = inertialNavigation(i,8:10)'*deg_to_rad; % roll,pitch,heading
    f_ib_b = inertialNavigation(i,11:end)'; 

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

    [r_caret_ea_e,v_caret_ea_e] = pv_NED_to_ECEF(gnssPosVel(i,2)*deg_to_rad,gnssPosVel(i,3)*deg_to_rad,gnssPosVel(i,4),gnssPosVel(i,5:7)');

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

    solution = [solution;[time,L_b*rad_to_deg,lambda_b*rad_to_deg,h_b,v_eb_n',eul']];

end

solutionTable = array2table(solution);
writetable(solutionTable,"task3Solution.csv","WriteVariableNames",false);