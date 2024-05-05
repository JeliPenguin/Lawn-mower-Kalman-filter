Define_Constants;
% Task 1a
phoneLat = -33.821075*deg_to_rad; %In radians
phoneLong = 151.188496*deg_to_rad; %In radians
phoneHeight = 120; %In metres
% [userPosPred,userVelPred] = pv_NED_to_ECEF(phoneLat,phoneLong,phoneHeight,0);
% disp(userPosPred)

% Task 1b
pseudoRanges = readtable("Workshop1_Pseudo_ranges.csv");
epochNum = height(pseudoRanges) - 1;
pseudoRangeRates = readtable("Workshop1_Pseudo_range_rates.csv");
loadedRow = table2array(pseudoRanges(1,:));
satelliteNumbers = loadedRow(2:end);

Omega_ie_e = [0,-omega_ie,0
              omega_ie,0,0
              0,0,0
              ];

userPosPred = [0;0;0]; % Initial user position prediction
for i=1:epochNum
    originalUserPosPred = userPosPred;
    % Load in data for current epoch
    loadedRow = table2array(pseudoRanges(i+1,:));
    time = loadedRow(1);
    measuredPseudoRanges = loadedRow(2:end);
    % disp(userPosPred)
    [res,userPosPred,outlierRemoval,maxJ] = processEpoch(satelliteNumbers,measuredPseudoRanges,userPosPred,time);
    % disp(outlierRemoval)
    if outlierRemoval
        measuredPseudoRanges(maxJ) = [];
        altSatelliteNumbers = satelliteNumbers;
        altSatelliteNumbers(maxJ) = [];
        [res,userPosPred,outlierRemoval,maxJ] = processEpoch(altSatelliteNumbers,measuredPseudoRanges,originalUserPosPred,time);
    end
    % disp(userPosPred)
    disp([res,outlierRemoval,maxJ])
end

function [res,userPosPred,outlierRemoval,maxJ] = processEpoch(satelliteNumbers,measuredPseudoRanges,userPosPred,time)
    Define_Constants;
    numberOfSatellites = length(satelliteNumbers);
    outlierRemoval = false;
    r__caret_es_e = zeros(3,numberOfSatellites);
    v_caret_ea_e = zeros(3,numberOfSatellites);
    v_caret_ea_e_minus = zeros(3,numberOfSatellites);
    r_caret_as_minus = zeros(1,numberOfSatellites);
    r_caret_dot_as_minus = zeros(1,numberOfSatellites);
    
    for i = 1:numberOfSatellites
        [sat_r_es_e,sat_v_es_e]= Satellite_position_and_velocity(time,satelliteNumbers(i));
        r__caret_es_e(:,i) = sat_r_es_e;
        v_caret_ea_e(:,i) = sat_v_es_e;

        diff = r__caret_es_e(:,i)-userPosPred;
        r_caret_as_minus(i) = sqrt(diff.'*diff); % Initial iteration assuming identity sagnac effect matrix
    end
    % disp(r_es_e)

    sigma_p = 5;
    T = 6; % Outlier detection threshold
    diff = 1;
    while diff > 0.0001
        r_caret_as_minus = calcRangePred(r__caret_es_e,userPosPred,r_caret_as_minus);
        H_e_G = calcMeasurementMatrix(r__caret_es_e,userPosPred,r_caret_as_minus);
        % disp(H_e_G)
        pred_receiver_clock_offset = 0;
        delta_z_min = zeros(numberOfSatellites,1); % Measurement Innovation vector
        for i=1:numberOfSatellites
            delta_z_min(i) = measuredPseudoRanges(i)-r_caret_as_minus(i)-pred_receiver_clock_offset;
        end

        x_caret_min = [userPosPred;pred_receiver_clock_offset];
        x_caret_plus = x_caret_min + inv(H_e_G.'*H_e_G)*H_e_G.'*delta_z_min;
        r_caret_e_plus_ea = x_caret_plus(1:3);
        % disp(r_caret_e_plus_ea)
        receiver_clock_offset = x_caret_plus(end);
        
        [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_caret_e_plus_ea,0);
        lat = L_b * rad_to_deg;
        long = lambda_b * rad_to_deg;
        % disp([lat,long,h_b])
        diff = norm(r_caret_e_plus_ea-userPosPred);
        % diff = diff + 1;
        userPosPred = r_caret_e_plus_ea; % For next iteration
    end

    % Outlier detection
    v = (H_e_G*inv(H_e_G.'*H_e_G)*H_e_G.'-eye(numberOfSatellites))*delta_z_min; % residuals vector
    C_v = (eye(numberOfSatellites)-H_e_G*inv(H_e_G.'*H_e_G)*H_e_G.')*sigma_p^2;
    maxJ = -1;
    maxResidual = -inf;
    for j=1:numberOfSatellites
        normalized_residual = norm(v(j))/sqrt(C_v(j,j));
        if normalized_residual > T
            % disp(["Outlier detected at: ",j," Residual: ",normalized_residual])
            % sprintf("Outlier detected with Satellite %d, Residual: %f",satelliteNumbers(j),normalized_residual)
            if normalized_residual > maxResidual
                maxJ = j;
                maxResidual = normalized_residual;
            end
            outlierRemoval = true;
        end
    end

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