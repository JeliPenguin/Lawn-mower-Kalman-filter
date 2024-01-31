% Load in data

pseudoRangeFile = readtable("Workshop2_Pseudo_ranges.csv");
epochNum = height(pseudoRangeFile) - 1;
pseudoRangeRateFile = readtable("Workshop2_Pseudo_range_rates.csv");
ecefPosFile = readtable("Workshop2_GNSS_Pos_ECEF.csv");

satelliteNumRow = table2array(pseudoRangeFile(1,:));
satelliteNumbers = satelliteNumRow(2:end);

pseudoRanges = table2array(pseudoRangeFile(2:end,2:end));
pseudoRangeRates = table2array(pseudoRangeRateFile(2:end,2:end));
ecefPos = table2array(ecefPosFile(1:end,2:end));

tau = 1; % propagation interval
S_a = 5; % acceleration power spectral density

x_caret_plus_k_m_1 = [2447019;-5884199;-284783;184;77;0];
P_plus_k_m_1 = eye(6) * 100;
P_plus_k_m_1(4,4) = 25;
P_plus_k_m_1(5,5) = 25;
P_plus_k_m_1(6,6) = 25;

for epoch=1:epochNum

    phi_k_m_1 = [eye(3),eye(3)*tau;zeros(3),eye(3)];
    
    Q_k_m_1 = [S_a*tau^3*eye(3)/3,S_a*tau^2*eye(3)/2;S_a*tau^2*eye(3)/2,S_a*tau*eye(3)];
    
    x_caret_minus_k = phi_k_m_1 * x_caret_plus_k_m_1;
    P_minus_k = phi_k_m_1*P_plus_k_m_1*phi_k_m_1.' + Q_k_m_1;
    
    H_k = [eye(3),zeros(3)];
    error_std = 2.5;
    R_k = eye(3) * error_std^2;
    
    % Kalman gain matrix
    K_k = P_minus_k * H_k.'*inv(H_k*P_minus_k*H_k.'+R_k);
    
    % measurement innovation
    r_tilde = ecefPos(epoch,:).';
    delta_z_min = r_tilde-x_caret_minus_k(1:3);
    
    x_caret_plus_k_m_1 = x_caret_minus_k + K_k * delta_z_min;
    P_plus_k_m_1 = (eye(6)-K_k*H_k)*P_minus_k;
    r_eb_e = x_caret_plus_k_m_1(1:3);
    v_eb_e = x_caret_plus_k_m_1(4:end);
    
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_eb_e,v_eb_e);
    L_b = L_b * 180 / pi;
    lambda_b = lambda_b * 180/pi;
    
    
    disp([epoch-1,L_b,lambda_b,h_b,v_eb_n.'])
end