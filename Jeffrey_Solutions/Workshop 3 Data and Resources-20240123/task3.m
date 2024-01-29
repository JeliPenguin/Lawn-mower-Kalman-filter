Define_Constants;

task2SolutionTable = readtable("task2Solution.csv");
gnssPosVel = table2array(task2SolutionTable);

inertiaNavigationTable = readtable("Aircraft_inertial.csv");
inertialNavigation = table2array(inertiaNavigationTable);
epochNum = height(inertialNavigation);

[x_caret_k_m_1_plus,P_plus_k_m_1] = Initialise_Integration_KF;

% Define Constants
tau_s = 1;

for i=1:1
    time = inertialNavigation(i,1);
    L_b = inertialNavigation(i,2) * deg_to_rad;
    lambda_b = inertialNavigation(i,3) * deg_to_rad;
    h_b = inertialNavigation(i,4);
    v_eb_n = inertialNavigation(i,5:7)';
    psi_nb = inertialNavigation(i,8:10)';
    f_ib_b = inertialNavigation(i,11:end)';

    delta_psi_eb_e = x_caret_k_m_1_plus(1:3);
    delta_v_eb_e = x_caret_k_m_1_plus(4:6);
    delta_r_eb_e = x_caret_k_m_1_plus(7:9);
    b_a = x_caret_k_m_1_plus(10:12);
    b_g = x_caret_k_m_1_plus(13:15);

    C_n_b = Euler_to_CTM(psi_nb);
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

    Q_k_m_1 = [LC_KF_config.gyro_noise_PSD*tau_s*eye(3),zeros(3),zeros(3),zeros(3),zeros(3)
               zeros(3),LC_KF_config.accel_noise_PSD*tau_s*eye(3),zeros(3),zeros(3),zeros(3)
               zeros(3),zeros(3),zeros(3),zeros(3),zeros(3)
               zeros(3),zeros(3),zeros(3),LC_KF_config.accel_bias_PSD*tau_s*eye(3),zeros(3)
               zeros(3),zeros(3),zeros(3),zeros(3),LC_KF_config.gyro_bias_PSD*tau_s*eye(3)
               ];

    x_caret_k_minus = Phi_k_m_1*x_caret_k_m_1_plus;
    P_k_minus = Phi_k_m_1*P_plus_k_m_1*Phi_k_m_1'+Q_k_m_1;


    H_k = [zeros(3),zeros(3),-eye(3),zeros(3),zeros(3)
           zeros(3),-eye(3),zeros(3),zeros(3),zeros(3)];

    R_k = [LC_KF_config.pos_meas_SD*eye(3),zeros(3)
           zeros(3),LC_KF_config.vel_meas_SD*eye(3)
           ];

    K_k = P_k_minus*H_k'*inv(H_k*P_k_minus*H_k'+R_k);


    [r_caret_ea_e,v_caret_ea_e,a] = NED_to_ECEF(gnssPosVel(i,2),gnssPosVel(i,3),gnssPosVel(i,4),gnssPosVel(i,5:7)',C_b_n);

    delta_z_k_minus = [r_caret_ea_e;v_caret_ea_e] - [r_eb_e;v_eb_e] + [delta_r_eb_e;delta_v_eb_e];

    x_caret_k_m_1_plus = x_caret_k_minus + K_k*delta_z_k_minus;
    P_plus_k_m_1 = (eye(15)-K_k*H_k)*P_k_minus;

    skewed = Skew_symmetric(delta_psi_eb_e);
    C_caret_b_e = (eye(3)-skewed)*C_b_e;
    v_caret_eb_e = v_eb_e - delta_v_eb_e;
    r_caret_eb_e = r_eb_e - delta_r_eb_e;

    [L_b,lambda_b,h_b,v_eb_n,C_b_n] = ECEF_to_NED(r_caret_eb_e,v_caret_eb_e,C_caret_b_e);
    eul = CTM_to_Euler(C_b_n);

    disp([L_b*rad_to_deg,lambda_b*rad_to_deg,h_b,v_eb_n',eul'*rad_to_deg])

end