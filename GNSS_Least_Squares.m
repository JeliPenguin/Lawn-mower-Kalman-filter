function [x_est,H_G_e,delta_z_min,delta_z_dot_min,r_caret_ea_e_minus,v_caret_ea_e_minus] = GNSS_Least_Squares(time,epoch_pseudo_range_measurements,epoch_pseudo_range_rate_measurements,satellite_numbers,r_caret_ea_e_minus,v_caret_ea_e_minus)
% Calculates least square solution for current epoch GNSS data
%
% Inputs:
%   epoch                                      Index of current epoch
%   times                                      Time of current epoch
%   epoch_pseudo_range_measurements            Array of pseudo-range measurements for all epochs from all satellites   
%   epoch_pseudo_range_rate_measurements       Array of pseudo-range rate measurements for all epochs from all satellites   
%   satellite_numbers                          Array of the number given to each of the satellites    
%   r_caret_ea_e_minus                         User position prediction
%   v_caret_ea_e_minus                         User velocity prediction
%
% Outputs:
%   x_est                                      Vector of initial state estimates:
%                                                Rows 1-3            estimated ECEF user position (m)
%                                                Rows 4-6            estimated ECEF user velocity (m/s)
%                                                Row 7               estimated receiver clock offset (m) 
%                                                Row 8               estimated receiver clock drift (m/s)
%  epoch_valid_indicies                        Array indicating whether the satellite
%                                              measurement is an outlier or not for current epoch
%  r_caret_ea_e_minus                          Position prediction for next epoch
%  v_caret_ea_e_minus                          Velocity prediction for next epoch

Define_Constants;

% Load in given data at current epoch
number_of_satellites = length(satellite_numbers);

% Initialization of range prediction array
r_caret_as_minus = zeros(1,number_of_satellites);

% Initialization of range rate prediction array
r_caret_dot_as_minus = zeros(1,number_of_satellites);

% Initialization of measurement Innovation vector for position
delta_z_min = zeros(number_of_satellites,1); 

% Initialization of measurement Innovation vector for velocity
delta_z_dot_min = zeros(number_of_satellites,1); 

r_caret_es_e = zeros(3,number_of_satellites);
v_caret_es_e = zeros(3,number_of_satellites);
H_G_e = zeros(number_of_satellites,4);

% Predicted receiver clock offset
delta_rho_caret_c_a_minus = 0; 

%Predicted receiver clock drift
delta_rho_caret_dot_c_a_minus = 0; 

for j = 1:number_of_satellites
    % Compute Cartesian ECEF positions and velocities of satellites at epoch 0
    [sat_r_caret_es_e,sat_v_es_e]= Satellite_position_and_velocity(time,satellite_numbers(j)); 
    r_caret_es_e(:,j) = sat_r_caret_es_e;
    v_caret_es_e(:,j) = sat_v_es_e;

    % Range Prediction initialization assuming identity sagnac effect matrix
    diff = r_caret_es_e(:,j)-r_caret_ea_e_minus;
    r_caret_as_minus(j) = sqrt(diff.'*diff);
end

diff = 1;

% Calculate least square solution for state estimate
while diff > 0.0001

    for j = 1:number_of_satellites
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
        delta_z_min(j) = epoch_pseudo_range_measurements(j)-r_caret_as_minus(j)-delta_rho_caret_c_a_minus;
        delta_z_dot_min(j) = epoch_pseudo_range_rate_measurements(j)-r_caret_dot_as_minus(j)-delta_rho_caret_dot_c_a_minus;
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

x_est = zeros(8,1);
x_est(1:3) = r_caret_ea_e_plus;
x_est(4:6) = v_caret_ea_e_plus;
x_est(7) = delta_rho_caret_c_a_plus;
x_est(8) = delta_rho_caret_dot_c_a_plus;

end  