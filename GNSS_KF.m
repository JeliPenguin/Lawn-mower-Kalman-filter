function gnss_solution = GNSS_KF(pseudo_range_data,pseudo_range_rate_data,times)
% Calculates GNSS Kalman Filter estimation, saving solutions in GNSS_Solution.csv
% Inputs:
%   pseudo_range_data                 Array of raw pseudo range measurement data
%   pseudo_range_rate_data            Array of raw pseudo range rate measurement data
%   times                             Array of time value for each epoch
%
% Outputs:
%   gnss_solution                     Calculated GNSS KF solution, for each row
%                                       Col 1: Time (s)
%                                       Col 2: Latitude (deg)
%                                       Col 3: Longitude (deg)
%                                       Col 4: Height (m)
%                                       Col 5: North Velocity (m/s)
%                                       Col 6: East Velocity (m/s)
%                                       Col 7: Down Velocity (m/s)

% Initializing constants
Define_Constants;

% Initialize GNSS KF Calculation constants
satellite_numbers = pseudo_range_data(1,2:end); % Array of satellite numbers
satellite_count = length(satellite_numbers); % Total count of satellites

epoch_num = height(times); % Total count of epochs

pseudo_ranges = pseudo_range_data(2:end,2:end);
pseudo_range_rates = pseudo_range_rate_data(2:end,2:end);

% Array for saving results
gnss_solution = [];

% Initialize initial state for Kalman Filter, giving initial x_est and P_matrix
[x_caret_plus_k_m_1,P_plus_k_m_1] = Initialise_GNSS_KF(times,pseudo_ranges,pseudo_range_rates,satellite_numbers,satellite_count,Omega_ie,sigma_p,T,sigma_co,sigma_cd,sigma_r,sigma_v);

% Compute transition matrix
phi_k_m_1 = CalcPhi(tau);

% Compute system noise covariance matrix
Q_k_m_1 = CalcQ(tau,S_a,S_cf,S_c_phi);

% Predicted psuedo ranges
r_caret_as_minus = zeros(satellite_count,1); 

for j = 1:satellite_count
    [sat_r_es_e,sat_v_es_e]= Satellite_position_and_velocity(times(1),satellite_numbers(j));
    diff = sat_r_es_e.'-x_caret_plus_k_m_1(1:3);
    r_caret_as_minus(j) = sqrt(diff.'*diff);
end % Predicted range initialization

u_e_a = zeros(3,satellite_count); % Line of sight
H_k = zeros(satellite_count*2,8);
H_k(1:satellite_count,7) = 1;
H_k(satellite_count+1:end,8) = 1;

r_caret_dot_as_minus = zeros(satellite_count,1); % Range rate predictions

delta_z_min = zeros(satellite_count*2,1); % Measurement innovation

for epoch = 1:epoch_num

    % Load data for current epoch
    time = times(epoch);
    current_epoch_pseudo_ranges = pseudo_ranges(epoch,:);
    current_epoch_pseudo_range_rates = pseudo_range_rates(epoch,:);

    % Propagate state estimate and error covariance
    x_caret_minus_k = phi_k_m_1 * x_caret_plus_k_m_1;
    P_minus_k = phi_k_m_1*P_plus_k_m_1*phi_k_m_1.' + Q_k_m_1;

    % Taking out each component of the state estimate vector
    r_caret_e_minus_ea_k = x_caret_minus_k(1:3);
    v_caret_e_minus_ea_k = x_caret_minus_k(4:6);
    roh_caret_a_minuc_c = x_caret_minus_k(7);
    roh_dot_caret_a_minuc_c = x_caret_minus_k(8);

    for j = 1:satellite_count
        C_I_e = CalcSagnacMatrix(r_caret_as_minus(j));

        % Range prediction from approximate user position to satellite
        [sat_r_es_e,sat_v_es_e]= Satellite_position_and_velocity(time,satellite_numbers(j));
        diff = C_I_e*sat_r_es_e.'-r_caret_e_minus_ea_k;
        r_caret_as_minus(j) = sqrt(diff.'*diff);

        % Line of sight vector calculation
        u_e_a(:,j) = diff / r_caret_as_minus(j);

        % Range rate prediction
        r_caret_dot_as_minus(j) = u_e_a(:,j).' * (C_I_e*(sat_v_es_e.'+Omega_ie*sat_r_es_e.')-(v_caret_e_minus_ea_k+Omega_ie*r_caret_e_minus_ea_k));
        
        % Measurement matrix calculation
        H_k(j,1:3) = -u_e_a(:,j);
        H_k(j+satellite_count,4:6) = -u_e_a(:,j);
    end
    
    % Measurement noise covariance matrix calculation
    R_k = CalcR(sigma_rho,sigma_rho_dot);
    
    % Kalman gain matrix calculation
    K_k = P_minus_k * H_k.'*inv(H_k*P_minus_k*H_k.'+R_k);
    
    % Measurement innovation calculation
    for j = 1:satellite_count
        delta_z_min(j) = current_epoch_pseudo_ranges(j) - r_caret_as_minus(j) - roh_caret_a_minuc_c;
        delta_z_min(j+satellite_count) = current_epoch_pseudo_range_rates(j) - r_caret_dot_as_minus(j) - roh_dot_caret_a_minuc_c;
    end

    % Update state estimates and error covariance
    x_caret_plus_k_m_1 = x_caret_minus_k + K_k * delta_z_min;
    P_plus_k_m_1 = (eye(8)-K_k*H_k)*P_minus_k;

    r_eb_e = x_caret_plus_k_m_1(1:3);
    v_eb_e = x_caret_plus_k_m_1(4:6);

    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_eb_e,v_eb_e);
    L_b = L_b * rad_to_deg;
    lambda_b = lambda_b * rad_to_deg;
    gnss_solution = [gnss_solution;time,L_b,lambda_b,h_b,v_eb_n.'];
end

outputTable = table(gnss_solution);
writetable(outputTable,"Solutions/GNSS_Solution"+".csv",'WriteVariableNames',0)
end