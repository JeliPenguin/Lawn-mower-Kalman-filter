function [x_est,P_matrix] = Initialise_GNSS_KF(times,pseudoRanges,pseudoRangeRates,satelliteNumbers,numberOfSatellites,Omega_ie,sigma_p,T,sigma_co,sigma_cd,sigma_r,sigma_v)
    %Initialise_GNSS_KF - Initializes the GNSS EKF state estimates and error
    %covariance matrix using one epoch least squares

    %   x_est                 Kalman filter estimates:
    %     Rows 1-3            estimated ECEF user position (m)
    %     Rows 4-6            estimated ECEF user velocity (m/s)
    %     Row 7               estimated receiver clock offset (m) 
    %     Row 8               estimated receiver clock drift (m/s)
    %   P_matrix              state estimation error covariance matrix
    
    % Load in given data at first epoch
    time = times(1);
    pseudoRangeMeasurements = pseudoRanges(1,:);
    pseudoRangeRateMeasurements = pseudoRangeRates(1,:);
    
    % Define variables
    r_caret_as_minus = zeros(1,numberOfSatellites);
    r_caret_dot_as_minus = zeros(1,numberOfSatellites);
    delta_z_min = zeros(numberOfSatellites,1); % Measurement Innovation vector
    delta_z_dot_min = zeros(numberOfSatellites,1); % Measurement Innovation vector
    r_caret_es_e = zeros(3,numberOfSatellites);
    v_caret_es_e = zeros(3,numberOfSatellites);
    H_G_e = zeros(numberOfSatellites,4);
    x_est = zeros(8,1);
    delta_rho_caret_c_a_minus = 0; % Predicted receiver clock offset
    delta_rho_caret_dot_c_a_minus = 0; %Predicted receiver clock drift
    
    r_caret_ea_e_minus = [0;0;0]; % Initial user position prediction
    v_caret_ea_e_minus = [0;0;0]; % Initial user velocity prediction
    
    for j = 1:numberOfSatellites
        [sat_r_caret_es_e,sat_v_es_e]= Satellite_position_and_velocity(time,satelliteNumbers(j));
        r_caret_es_e(:,j) = sat_r_caret_es_e;
        v_caret_es_e(:,j) = sat_v_es_e;
    
        % Position Prediction initialization assuming identity sagnac effect matrix
        diff = r_caret_es_e(:,j)-r_caret_ea_e_minus;
        r_caret_as_minus(j) = sqrt(diff.'*diff);
    end
    
    diff = 1;
    
    % Calculate least square solution for state estimate
    while diff > 0.0001
    
        for j = 1:numberOfSatellites
            C_I_e = CalcSagnacMatrix(r_caret_as_minus(j));
            diff = C_I_e*r_caret_es_e(:,j)-r_caret_ea_e_minus;
            % Actual calculation for the ranges prediction
            r_caret_as_minus(j) = sqrt(diff.'*diff);
    
            u_caret_aj_e_minus = (C_I_e*r_caret_es_e(:,j)-r_caret_ea_e_minus)/r_caret_as_minus(j);
            a = C_I_e*(v_caret_es_e(:,j)+Omega_ie*r_caret_es_e(:,j));
            b = v_caret_ea_e_minus + Omega_ie*r_caret_ea_e_minus;
            r_caret_dot_as_minus(j) = u_caret_aj_e_minus' * (a-b);
    
            % Calculate Measurement Matrix
            H_G_e(j,1:3) = -u_caret_aj_e_minus.';
            H_G_e(j,4) = 1;
    
            % Calculate innovation vectors
            delta_z_min(j) = pseudoRangeMeasurements(j)-r_caret_as_minus(j)-delta_rho_caret_c_a_minus;
    
            delta_z_dot_min(j) = pseudoRangeRateMeasurements(j)-r_caret_dot_as_minus(j)-delta_rho_caret_dot_c_a_minus;
        end 
    
        x_caret_min = [r_caret_ea_e_minus;delta_rho_caret_c_a_minus];
        x_caret_plus = x_caret_min + inv(H_G_e.'*H_G_e)*H_G_e.'*delta_z_min;
        r_caret_ea_e_plus = x_caret_plus(1:3);
        delta_rho_caret_c_a_plus = x_caret_plus(end);
    
        x_caret_min = [v_caret_ea_e_minus;delta_rho_caret_dot_c_a_minus];
        x_caret_plus = x_caret_min + inv(H_G_e.'*H_G_e)*H_G_e.'*delta_z_dot_min;
        v_caret_ea_e_plus = x_caret_plus(1:3);
        delta_rho_caret_dot_c_a_plus = x_caret_plus(end);
        
        diff = norm([r_caret_ea_e_plus,v_caret_ea_e_plus]-[r_caret_ea_e_minus,v_caret_ea_e_minus]);
    
        % For next iteration
        r_caret_ea_e_minus = r_caret_ea_e_plus;
        v_caret_ea_e_minus = v_caret_ea_e_plus;
    end
    
    hasOutlier = OutlierDetection(numberOfSatellites,H_G_e,delta_z_min,sigma_p,T);
    if hasOutlier
        disp("Outlier Detected")
    else
        disp("No Outliers Detected")
    end
    
    x_est(1:3) = r_caret_ea_e_plus;
    x_est(4:6) = v_caret_ea_e_plus;
    x_est(7) = delta_rho_caret_c_a_plus;
    x_est(8) = delta_rho_caret_dot_c_a_plus;
    
    % Initialise error covariance matrix
    P_matrix =  zeros(8);
    P_matrix(1,1) = sigma_r^2;
    P_matrix(2,2) = sigma_r^2;
    P_matrix(3,3) = sigma_r^2;
    P_matrix(4,4) = sigma_v^2;
    P_matrix(5,5) = sigma_v^2;
    P_matrix(6,6) = sigma_v^2;
    P_matrix(7,7) = sigma_co^2;
    P_matrix(8,8) = sigma_cd^2;
    
    % Ends
end