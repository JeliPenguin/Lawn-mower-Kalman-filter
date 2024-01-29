Define_Constants;
% Load in data

pseudoRangeFile = readtable("Workshop2_Pseudo_ranges.csv");
epochNum = height(pseudoRangeFile) - 1;
pseudoRangeRateFile = readtable("Workshop2_Pseudo_range_rates.csv");
ecefPosFile = readtable("Workshop2_GNSS_Pos_ECEF.csv");

satelliteNumRow = table2array(pseudoRangeFile(1,:));
satelliteNumbers = satelliteNumRow(2:end);
numOfSatellites = width(satelliteNumbers);

pseudoRanges = table2array(pseudoRangeFile(2:end,2:end));
pseudoRangeRates = table2array(pseudoRangeRateFile(2:end,2:end));
ecefPos = table2array(ecefPosFile(1:end,2:end));

[x_caret_plus_k_m_1,P_plus_k_m_1] = Initialise_GNSS_KF;

% disp(x_caret_plus_k_m_1)

r_caret_minus = zeros(numOfSatellites,1); % Predicted ranges
for j = 1:numOfSatellites
    [sat_r_es_e,sat_v_es_e]= Satellite_position_and_velocity(0,satelliteNumbers(j));
    diff = sat_r_es_e.'-x_caret_plus_k_m_1(1:3);
    r_caret_minus(j) = sqrt(diff.'*diff);
end % Predicted range initialization

u_e_a = zeros(3,numOfSatellites); % Line of sight
H_k = zeros(numOfSatellites*2,8);
H_k(1:numOfSatellites,7) = 1;
H_k(numOfSatellites+1:end,8) = 1;

rangeRatePred = zeros(numOfSatellites,1); % Range rate predictions

delta_z_min = zeros(numOfSatellites*2,1); % Measurement innovation

Define_Constants;
tau = 1; % propagation interval
S_a = 5; % acceleration power spectral density
S_c_phi = 0.01; % clock phase PSD
S_cf = 0.04; % clock frequency PSD

phi_k_m_1 = [eye(3),eye(3)*tau,zeros(3,1),zeros(3,1)
             zeros(3),eye(3),zeros(3,1),zeros(3,1)
             zeros(1,3),zeros(1,3),1,tau
             zeros(1,3),zeros(1,3),0,1];

Q_k_m_1 = [S_a*tau^3*eye(3)/3,S_a*tau^2*eye(3)/2,zeros(3,1),zeros(3,1)
           S_a*tau^2*eye(3)/2,S_a*tau*eye(3),zeros(3,1),zeros(3,1)
           zeros(1,3),zeros(1,3),S_c_phi*tau+S_cf*tau^3/3,S_cf*tau^2/2
           zeros(1,3),zeros(1,3),S_cf*tau^2/2,S_cf*tau
           ];

% disp(Q_k_m_1)
for epoch = 1:epochNum
% for epoch =1:1
    currentEpochPseudoRanges = pseudoRanges(epoch,:);
    currentEpochPseudoRangeRates = pseudoRangeRates(epoch,:);
    x_caret_minus_k = phi_k_m_1 * x_caret_plus_k_m_1;
    r_caret_e_minus_ea_k = x_caret_minus_k(1:3);
    v_caret_e_minus_ea_k = x_caret_minus_k(4:6);
    roh_caret_a_minuc_c = x_caret_minus_k(7);
    roh_dot_caret_a_minuc_c = x_caret_minus_k(8);
    P_minus_k = phi_k_m_1*P_plus_k_m_1*phi_k_m_1.' + Q_k_m_1;
    % disp(x_caret_minus_k) % Part D
    % disp(P_minus_k) % Part E
    for j = 1:numOfSatellites
        C_I_e = eye(3);
        C_I_e(1,2) = omega_ie*r_caret_minus(j)/c;
        C_I_e(2,1) = -omega_ie*r_caret_minus(j)/c;
        [sat_r_es_e,sat_v_es_e]= Satellite_position_and_velocity(epoch-1,satelliteNumbers(j));
        diff = C_I_e*sat_r_es_e.'-r_caret_e_minus_ea_k;
        r_caret_minus(j) = sqrt(diff.'*diff);
        u_e_a(:,j) = diff / r_caret_minus(j);
        rangeRatePred(j) = u_e_a(:,j).' * (C_I_e*(sat_v_es_e.'+Omega_ie*sat_r_es_e.')-(v_caret_e_minus_ea_k+Omega_ie*r_caret_e_minus_ea_k));
        H_k(j,1:3) = -u_e_a(:,j);
        H_k(j+numOfSatellites,4:6) = -u_e_a(:,j);
    end

    sigma_p = 10;
    sigma_r = 0.05;
    R_k = [eye(10,10) * sigma_p^2,zeros(10,10);zeros(10,10),eye(10,10)*sigma_r^2];
    
    % Kalman gain matrix
    K_k = P_minus_k * H_k.'*inv(H_k*P_minus_k*H_k.'+R_k);
    % disp(K_k) % Part k
    
    % measurement innovation
    for j = 1:numOfSatellites
        delta_z_min(j) = currentEpochPseudoRanges(j) - r_caret_minus(j) - roh_caret_a_minuc_c;
        delta_z_min(j+numOfSatellites) = currentEpochPseudoRangeRates(j) - rangeRatePred(j) - roh_dot_caret_a_minuc_c;
    end

    % disp(delta_z_min)

    x_caret_plus_k_m_1 = x_caret_minus_k + K_k * delta_z_min;
    % disp(x_caret_plus_k_m_1)
    P_plus_k_m_1 = (eye(8)-K_k*H_k)*P_minus_k;
    % disp(P_plus_k_m_1)
    r_eb_e = x_caret_plus_k_m_1(1:3);
    v_eb_e = x_caret_plus_k_m_1(4:6);

    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_eb_e,v_eb_e);
    L_b = L_b * rad_to_deg;
    lambda_b = lambda_b * rad_to_deg;
    disp([epoch-1,L_b,lambda_b,h_b,v_eb_n.'])
    % disp(epoch)
end