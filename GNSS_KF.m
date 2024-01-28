Define_Constants;
Load_Data;

GNSS_Solution = [];

% Initialize initial state for Kalman Filter, giving initial x_est and
% P_matrix
Initialise_GNSS_KF;

x_caret_plus_k_m_1 = x_est;
P_plus_k_m_1 = P_matrix;

% Compute transition matrix
phi_k_m_1 = [eye(3),eye(3)*tau,zeros(3,1),zeros(3,1)
             zeros(3),eye(3),zeros(3,1),zeros(3,1)
             zeros(1,3),zeros(1,3),1,tau
             zeros(1,3),zeros(1,3),0,1];

% Compute system noise covariance matrix
Q_k_m_1 = [S_a*tau^3*eye(3)/3,S_a*tau^2*eye(3)/2,zeros(3,1),zeros(3,1)
           S_a*tau^2*eye(3)/2,S_a*tau*eye(3),zeros(3,1),zeros(3,1)
           zeros(1,3),zeros(1,3),S_c_phi*tau+S_cf*tau^3/3,S_cf*tau^2/2
           zeros(1,3),zeros(1,3),S_cf*tau^2/2,S_cf*tau
           ];

r_caret_as_minus = zeros(numberOfSatellites,1); % Predicted psuedo ranges
for j = 1:numberOfSatellites
    [sat_r_es_e,sat_v_es_e]= Satellite_position_and_velocity(times(1),satelliteNumbers(j));
    diff = sat_r_es_e.'-x_caret_plus_k_m_1(1:3);
    r_caret_as_minus(j) = sqrt(diff.'*diff);
end % Predicted range initialization

u_e_a = zeros(3,numberOfSatellites); % Line of sight
H_k = zeros(numberOfSatellites*2,8);
H_k(1:numberOfSatellites,7) = 1;
H_k(numberOfSatellites+1:end,8) = 1;

r_caret_dot_as_minus = zeros(numberOfSatellites,1); % Range rate predictions

delta_z_min = zeros(numberOfSatellites*2,1); % Measurement innovation

for epoch = 1:epochNum

    % Load data for current epoch
    time = times(epoch);
    currentEpochPseudoRanges = pseudoRanges(epoch,:);
    currentEpochPseudoRangeRates = pseudoRangeRates(epoch,:);

    % Propagate state estimate and error covariance
    x_caret_minus_k = phi_k_m_1 * x_caret_plus_k_m_1;
    P_minus_k = phi_k_m_1*P_plus_k_m_1*phi_k_m_1.' + Q_k_m_1;

    % Taking out each component of the state estimate vector
    r_caret_e_minus_ea_k = x_caret_minus_k(1:3);
    v_caret_e_minus_ea_k = x_caret_minus_k(4:6);
    roh_caret_a_minuc_c = x_caret_minus_k(7);
    roh_dot_caret_a_minuc_c = x_caret_minus_k(8);

    for j = 1:numberOfSatellites
        C_I_e = CalcSagnacMatrix(r_caret_as_minus(j));

        % Range prediction from approximate user position to satellite
        [sat_r_es_e,sat_v_es_e]= Satellite_position_and_velocity(time,satelliteNumbers(j));
        diff = C_I_e*sat_r_es_e.'-r_caret_e_minus_ea_k;
        r_caret_as_minus(j) = sqrt(diff.'*diff);

        % Line of sight vector calculation
        u_e_a(:,j) = diff / r_caret_as_minus(j);

        % Range rate prediction
        r_caret_dot_as_minus(j) = u_e_a(:,j).' * (C_I_e*(sat_v_es_e.'+Omega_ie*sat_r_es_e.')-(v_caret_e_minus_ea_k+Omega_ie*r_caret_e_minus_ea_k));
        
        % Measurement matrix calculation
        H_k(j,1:3) = -u_e_a(:,j);
        H_k(j+numberOfSatellites,4:6) = -u_e_a(:,j);
    end
    
    % Measurement noise covariance matrix calculation
    R_k = [eye(10,10) * sigma_rho^2,zeros(10,10);zeros(10,10),eye(10,10)*sigma_rho_dot^2];
    
    % Kalman gain matrix calculation
    K_k = P_minus_k * H_k.'*inv(H_k*P_minus_k*H_k.'+R_k);
    
    % Measurement innovation calculation
    for j = 1:numberOfSatellites
        delta_z_min(j) = currentEpochPseudoRanges(j) - r_caret_as_minus(j) - roh_caret_a_minuc_c;
        delta_z_min(j+numberOfSatellites) = currentEpochPseudoRangeRates(j) - r_caret_dot_as_minus(j) - roh_dot_caret_a_minuc_c;
    end

    % Update state estimates and error covariance
    x_caret_plus_k_m_1 = x_caret_minus_k + K_k * delta_z_min;
    P_plus_k_m_1 = (eye(8)-K_k*H_k)*P_minus_k;

    r_eb_e = x_caret_plus_k_m_1(1:3);
    v_eb_e = x_caret_plus_k_m_1(4:6);

    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_eb_e,v_eb_e);
    L_b = L_b * rad_to_deg;
    lambda_b = lambda_b * rad_to_deg;
    GNSS_Solution = [GNSS_Solution;time,L_b,lambda_b,h_b,v_eb_n.'];
end

outputTable = table(GNSS_Solution);
writetable(outputTable,"GNSS_Solution"+".csv",'WriteVariableNames',0)