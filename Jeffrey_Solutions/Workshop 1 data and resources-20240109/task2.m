Define_Constants;
% Task 1a
% phoneLat = -33.821075*deg_to_rad; %In radians
% phoneLong = 151.188496*deg_to_rad; %In radians
% phoneHeight = 120; %In metres
% % [userPosPred,userVelPred] = pv_NED_to_ECEF(phoneLat,phoneLong,phoneHeight,0);
% % disp(userPosPred)

% Task 1b
pseudoRanges = readtable("Workshop1_Pseudo_ranges.csv");
epochNum = height(pseudoRanges) - 1;
pseudoRangeRates = readtable("Workshop1_Pseudo_range_rates.csv");
loadedRow = table2array(pseudoRanges(1,:));
satelliteNumbers = loadedRow(2:end);

r_caret_ea_e_minus = [0;0;0]; % Initial user position prediction
v_caret_ea_e_minus = [0;0;0];
for i=1:epochNum
    originalUserPosPred = r_caret_ea_e_minus;
    % Load in data for current epoch
    loadedRow = table2array(pseudoRanges(i+1,:));
    loadedRowRanges = table2array(pseudoRangeRates(i+1,:));
    time = loadedRow(1);
    measuredPseudoRanges = loadedRow(2:end);
    measuredPseudoRangeRates = loadedRowRanges(2:end);
    % disp(userPosPred)
    [res,r_caret_ea_e_minus] = processEpoch(satelliteNumbers,measuredPseudoRanges,r_caret_ea_e_minus,time,measuredPseudoRangeRates,v_caret_ea_e_minus);
end

function [res,r_caret_ea_e_minus] = processEpoch(satelliteNumbers,measuredPseudoRanges,r_caret_ea_e_minus,time,measuredPseudoRangeRates,v_caret_ea_e_minus)
    Define_Constants;
    numberOfSatellites = length(satelliteNumbers);
    r_caret_es_e = zeros(3,numberOfSatellites);
    v_caret_ea_e = zeros(3,numberOfSatellites);
    
    r_caret_as_minus = zeros(1,numberOfSatellites);
    r_caret_dot_as_minus = zeros(1,numberOfSatellites);
    
    for i = 1:numberOfSatellites
        [sat_r_es_e,sat_v_es_e]= Satellite_position_and_velocity(time,satelliteNumbers(i));
        r_caret_es_e(:,i) = sat_r_es_e;
        v_caret_ea_e(:,i) = sat_v_es_e;
        diff = r_caret_es_e(:,i)-r_caret_ea_e_minus;
        r_caret_as_minus(i) = sqrt(diff.'*diff); % Initial iteration assuming identity sagnac effect matrix
    end
    % disp(r_es_e)
    sigma_p = 5;
    T = 6; % Outlier detection threshold
    diff = 1;
    while diff > 0.00001
        for i =1:numberOfSatellites
            C_I_e = calcSagnacMatrix(r_caret_as_minus(i));
            diff = C_I_e*r_caret_es_e(:,i)-r_caret_ea_e_minus;
            r_caret_as_minus(i) = sqrt(diff.'*diff);
            
            u_caret_aj_e_minus = (C_I_e*r_caret_es_e(:,i)-r_caret_ea_e_minus)/r_caret_as_minus(i);
            a = C_I_e*(v_caret_ea_e(:,i)+Omega_ie*r_caret_es_e(:,i));
            b = v_caret_ea_e_minus + Omega_ie*r_caret_ea_e_minus;
            r_caret_dot_as_minus(i) = u_caret_aj_e_minus' * (a-b);
        end
        H_e_G = calcMeasurementMatrix(r_caret_es_e,r_caret_ea_e_minus,r_caret_as_minus);
        % disp(H_e_G)
        pred_receiver_clock_offset = 0;
        delta_z_min = zeros(numberOfSatellites,1); % Measurement Innovation vector
        for i=1:numberOfSatellites
            delta_z_min(i) = measuredPseudoRanges(i)-r_caret_as_minus(i)-pred_receiver_clock_offset;
        end

        x_caret_min = [r_caret_ea_e_minus;pred_receiver_clock_offset];
        x_caret_plus = x_caret_min + inv(H_e_G.'*H_e_G)*H_e_G.'*delta_z_min;
        r_caret_e_plus_ea = x_caret_plus(1:3);
        % disp(r_caret_e_plus_ea)
        receiver_clock_offset = x_caret_plus(end);

        for i=1:numberOfSatellites
            delta_z_min(i) = measuredPseudoRangeRates(i)-r_caret_dot_as_minus(i)-pred_receiver_clock_offset;
        end

        x_caret_min = [v_caret_ea_e_minus;pred_receiver_clock_offset];
        x_caret_plus = x_caret_min + inv(H_e_G.'*H_e_G)*H_e_G.'*delta_z_min;
        v_caret_ea_e_plus = x_caret_plus(1:3);

        
        [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_caret_e_plus_ea,v_caret_ea_e_plus);
        lat = L_b * rad_to_deg;
        long = lambda_b * rad_to_deg;
        diff = norm(r_caret_e_plus_ea-r_caret_ea_e_minus) + norm(v_caret_ea_e_plus-v_caret_ea_e_minus);
        r_caret_ea_e_minus = r_caret_e_plus_ea; % For next iteration
        v_caret_ea_e_minus = v_caret_ea_e_plus;
    end
    disp([time,v_eb_n'])
    res = [time,lat,long,h_b];
end

function H_e_G = calcMeasurementMatrix(satPositions,userPosPred,r_caret_min_aj) 
    numberOfSatellites = length(satPositions);
    H_e_G = zeros(numberOfSatellites,4);
    for i=1:numberOfSatellites
        C_I_e = calcSagnacMatrix(r_caret_min_aj(i));
        line_of_sight_vec = (C_I_e*satPositions(:,i)-userPosPred)/r_caret_min_aj(i);
        H_e_G(i,1:3) = -line_of_sight_vec.';
        H_e_G(i,4) = 1;
    end
end

function ret_r_caret_min_aj = calcRangePred(satPositions,userPosPred,r_caret_min_aj)
    numberOfSatellites = length(satPositions);
    ret_r_caret_min_aj = zeros(1,numberOfSatellites);
    for i = 1:numberOfSatellites
        C_I_e = calcSagnacMatrix(r_caret_min_aj(i));
        diff = C_I_e*satPositions(:,i)-userPosPred;
        ret_r_caret_min_aj(i) = sqrt(diff.'*diff);
    end % Actual calculation for the ranges prediction
end

function C_I_e = calcSagnacMatrix(range)
    Define_Constants;
    rep = omega_ie * range/c;
    C_I_e = [1,rep,0;-rep,1,0;0,0,1];
end