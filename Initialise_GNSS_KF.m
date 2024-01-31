function [x_est,P_matrix] = Initialise_GNSS_KF(times,pseudoRanges,pseudoRangeRates,satelliteNumbers,numberOfSatellites)
%Initialise_GNSS_KF - Initializes the GNSS EKF state estimates and error
%covariance matrix using one epoch least squares
%
% Inputs:
%   times                    Array of time value for each epoch
%   pseudoRanges             Array of pseudo-range measurements for all epochs from all satellites   
%   pseudoRangeRates         Array of pseudo-range rate measurements for all epochs from all satellites   
%   satelliteNumbers         Array of the number given to each of the satellites       
%   numberOfSatellites       Total count of number of satellites giving the measurements                    
%
% Outputs:
%   x_est                    Vector of initial state estimates:
%                               Rows 1-3            estimated ECEF user position (m)
%                               Rows 4-6            estimated ECEF user velocity (m/s)
%                               Row 7               estimated receiver clock offset (m) 
%                               Row 8               estimated receiver clock drift (m/s)
%   P_matrix                 State estimation error covariance matrix

Define_Constants;

% Load in given data at first epoch
pseudoRangeMeasurements = pseudoRanges(1,:);
pseudoRangeRateMeasurements = pseudoRangeRates(1,:);

% Initialization of range prediction array
r_caret_as_minus = zeros(1,numberOfSatellites);

% Initialization of range rate prediction array
r_caret_dot_as_minus = zeros(1,numberOfSatellites);

% Initialization of measurement Innovation vector for position
delta_z_min = zeros(numberOfSatellites,1); 

% Initialization of measurement Innovation vector for velocity
delta_z_dot_min = zeros(numberOfSatellites,1); 

r_caret_es_e = zeros(3,numberOfSatellites);
v_caret_es_e = zeros(3,numberOfSatellites);
H_G_e = zeros(numberOfSatellites,4);

% Predicted receiver clock offset
delta_rho_caret_c_a_minus = 0; 

%Predicted receiver clock drift
delta_rho_caret_dot_c_a_minus = 0; 

% Initial user position prediction
r_caret_ea_e_minus = [0;0;0]; 

% Initial user velocity prediction
v_caret_ea_e_minus = [0;0;0]; 

for j = 1:numberOfSatellites
    % Compute Cartesian ECEF positions and velocities of satellites at epoch 0
    [sat_r_caret_es_e,sat_v_es_e]= Satellite_position_and_velocity(times(1),satelliteNumbers(j)); 
    r_caret_es_e(:,j) = sat_r_caret_es_e;
    v_caret_es_e(:,j) = sat_v_es_e;

    % Range Prediction initialization assuming identity sagnac effect matrix
    diff = r_caret_es_e(:,j)-r_caret_ea_e_minus;
    r_caret_as_minus(j) = sqrt(diff.'*diff);
end

diff = 1;

% Calculate least square solution for state estimate
while diff > 0.0001

    for j = 1:numberOfSatellites
        % Caculate Sagnac matrix given range prediction of current satellite
        C_I_e = CalcSagnacMatrix(r_caret_as_minus(j));
        diff = C_I_e*r_caret_es_e(:,j)-r_caret_ea_e_minus;
        
        % Calculation for the ranges prediction of current satellite
        r_caret_as_minus(j) = sqrt(diff.'*diff);

        % Compute line-of-sight unit vector
        u_caret_aj_e_minus = (C_I_e*r_caret_es_e(:,j)-r_caret_ea_e_minus)/r_caret_as_minus(j);

        % Calculation for the range rate prediction of current satellite
        a = C_I_e*(v_caret_es_e(:,j)+Omega_ie*r_caret_es_e(:,j));
        b = v_caret_ea_e_minus + Omega_ie*r_caret_ea_e_minus;
        r_caret_dot_as_minus(j) = u_caret_aj_e_minus' * (a-b);

        % Calculating Measurement Matrix
        H_G_e(j,1:3) = -u_caret_aj_e_minus.';
        H_G_e(j,4) = 1;

        % Calculate innovation vectors for position and velocity
        delta_z_min(j) = pseudoRangeMeasurements(j)-r_caret_as_minus(j)-delta_rho_caret_c_a_minus;
        delta_z_dot_min(j) = pseudoRangeRateMeasurements(j)-r_caret_dot_as_minus(j)-delta_rho_caret_dot_c_a_minus;
    end 

    % Formulate predicted state vector for position
    x_caret_min = [r_caret_ea_e_minus;delta_rho_caret_c_a_minus];

    % Unweighted least-squares calculation, computing position and receiver
    % clock offset
    x_caret_plus = x_caret_min + inv(H_G_e.'*H_G_e)*H_G_e.'*delta_z_min;
    r_caret_ea_e_plus = x_caret_plus(1:3);
    delta_rho_caret_c_a_plus = x_caret_plus(end);
    
    % Formulate predicted state vector for velocity
    x_caret_min = [v_caret_ea_e_minus;delta_rho_caret_dot_c_a_minus];
    % Unweighted least-squares calculation, computing velocity and receiver
    % clock drift
    x_caret_plus = x_caret_min + inv(H_G_e.'*H_G_e)*H_G_e.'*delta_z_dot_min;
    v_caret_ea_e_plus = x_caret_plus(1:3);
    delta_rho_caret_dot_c_a_plus = x_caret_plus(end);

    % Calculate difference compared to previous computation
    diff = norm([r_caret_ea_e_plus,v_caret_ea_e_plus]-[r_caret_ea_e_minus,v_caret_ea_e_minus]);

    % For next iteration
    r_caret_ea_e_minus = r_caret_ea_e_plus;
    v_caret_ea_e_minus = v_caret_ea_e_plus;
end

% Outlier Detection
hasOutlier = OutlierDetection(numberOfSatellites,H_G_e,delta_z_min,sigma_p,T);
if hasOutlier
    disp("Outlier Detected")
else
    disp("No Outliers Detected")
end

% Setting return values for x_est and P_matrix
x_est = zeros(8,1);
x_est(1:3) = r_caret_ea_e_plus;
x_est(4:6) = v_caret_ea_e_plus;
x_est(7) = delta_rho_caret_c_a_plus;
x_est(8) = delta_rho_caret_dot_c_a_plus;

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