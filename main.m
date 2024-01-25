Define_Constants;
Load_Data;

% Configure outputs
outputFile = "Output_Profile";
% 1: Time
% 2: Geodetic latitude in degrees
% 3: Geodetic longitude in degrees 
% 4: North velocity in m/s
% 5: East velocity in m/s
% 6: Heading in degrees.
outputs = zeros(epochNum,6);

% Define additional constants
S_c_phi = 0.01; % clock phase PSD
S_cf = 0.04; % clock frequency PSD
sigma_p = 5;
T = 6; % Outlier detection threshold

% Define variables
r_caret_as_minus = zeros(1,numberOfSatellites);
r_caret_dot_as_minus = zeros(1,numberOfSatellites);
delta_z_min = zeros(numberOfSatellites,1); % Measurement Innovation vector
delta_z_dot_min = zeros(numberOfSatellites,1); % Measurement Innovation vector
r_caret_es_e = zeros(3,numberOfSatellites);
v_caret_es_e = zeros(3,numberOfSatellites);
H_G_e = zeros(numberOfSatellites,4);

r_caret_ea_e_minus = [0;0;0]; % Initial user position prediction
v_caret_ea_e_minus = [0;0;0]; % Initial user velocity prediction

delta_rho_caret_c_a_minus = 0; % Predicted receiver clock offset
delta_rho_caret_dot_c_a_minus = 0; %Predicted receiver clock drift

for i=1:epochNum
% for i=1:5
    % Load data for current epoch
    time = times(i);
    pseudoRangeMeasurements = pseudoRanges(i,:);
    pseudoRangeRateMeasurements = pseudoRangeRates(i,:);

    for j = 1:numberOfSatellites
        [sat_r_caret_es_e,sat_v_es_e]= Satellite_position_and_velocity(time,satelliteNumbers(j));
        r_caret_es_e(:,j) = sat_r_caret_es_e;
        v_caret_es_e(:,j) = sat_v_es_e;
    
        % Position Prediction initialization assuming identity sagnac effect matrix
        diff = r_caret_es_e(:,j)-r_caret_ea_e_minus;
        r_caret_as_minus(j) = sqrt(diff.'*diff);
    end

    diff = 1;
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
        
        [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_caret_ea_e_plus,v_caret_ea_e_plus);
        lat = L_b * rad_to_deg;
        long = lambda_b * rad_to_deg;
        diff = norm(r_caret_ea_e_plus-r_caret_ea_e_minus) + norm(v_caret_ea_e_plus-v_caret_ea_e_minus);

        % For next iteration
        r_caret_ea_e_minus = r_caret_ea_e_plus;
        v_caret_ea_e_minus = v_caret_ea_e_plus;
    end
    % disp([time,lat,long,h_b])

    heading = 0; % in radians
    outputs(i,:) = [time,lat,long,v_eb_n(1),v_eb_n(2),heading];
end

% Write to table and export
outputTable = table(outputs);
writetable(outputTable,outputFile+".csv",'WriteVariableNames',0)
