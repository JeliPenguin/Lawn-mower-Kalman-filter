function [x_est,P_matrix,valid_pseudo_range_indicies,valid_pseudo_range_rate_indicies] = Initialise_GNSS_KF(times,pseudo_ranges,pseudo_range_rates,satellite_numbers)
%Initialise_GNSS_KF - Initializes the GNSS EKF state estimates and error
%covariance matrix using one epoch least squares
%
% Inputs:
%   times                               Array of time value for each epoch
%   pseudo_ranges                       Array of pseudo-range measurements for all epochs from all satellites   
%   pseudo_range_rates                  Array of pseudo-range rate measurements for all epochs from all satellites   
%   satellite_numbers                   Array of the number given to each of the satellites       
%   number_of_satellites                Total count of number of satellites giving the measurements                    
%
% Outputs:
%   x_est                               Vector of initial state estimates:
%                                           Rows 1-3            estimated ECEF user position (m)
%                                           Rows 4-6            estimated ECEF user velocity (m/s)
%                                           Row 7               estimated receiver clock offset (m) 
%                                           Row 8               estimated receiver clock drift (m/s)
%   P_matrix                            State estimation error covariance matrix
%   valid_pseudo_range_indicies         Matrix indicating whether the satellite pseudo range
%                                       is an outlier or not (number of epochs X number of satellites)
%   valid_pseudo_range_rate_indicies    Matrix indicating whether the satellite pseudo range rate
%                                       is an outlier or not (number of epochs X number of satellites)

Define_Constants;

% Initial user position prediction
r_caret_ea_e_minus = [0;0;0]; 

% Initial user velocity prediction
v_caret_ea_e_minus = [0;0;0]; 

epoch_num = height(times);

valid_pseudo_range_indicies = ones(epoch_num,length(satellite_numbers));
valid_pseudo_range_rate_indicies = ones(epoch_num,length(satellite_numbers));

% First epoch outlier correction

original_r_caret_ea_e_minus = r_caret_ea_e_minus;
original_v_caret_ea_e_minus = v_caret_ea_e_minus;
valid_satellite_numbers = satellite_numbers;
epoch_pseudo_range_measurements = pseudo_ranges(1,:);
epoch_pseudo_range_rate_measurements = pseudo_range_rates(1,:);
[x_est_m,H_G_e,delta_z_min,delta_z_dot_min,r_caret_ea_e_minus,v_caret_ea_e_minus] = GNSS_Least_Squares(times(1),epoch_pseudo_range_measurements,epoch_pseudo_range_rate_measurements,satellite_numbers,r_caret_ea_e_minus,v_caret_ea_e_minus);

[~,maxJ] = Outlier_Detection(length(satellite_numbers),H_G_e,delta_z_min,sigma_rho);
[~,maxJDot] = Outlier_Detection(length(satellite_numbers),H_G_e,delta_z_dot_min,sigma_rho_dot);

while maxJ > 0
    % Do the calculation again with outlier removed
    epoch_pseudo_range_measurements(maxJ) = [];
    epoch_pseudo_range_rate_measurements(maxJ) = [];
    valid_satellite_numbers(maxJ) = [];
    [x_est_m,H_G_e,delta_z_min,delta_z_dot_min,r_caret_ea_e_minus,v_caret_ea_e_minus] = GNSS_Least_Squares(times(1),epoch_pseudo_range_measurements,epoch_pseudo_range_rate_measurements,valid_satellite_numbers,original_r_caret_ea_e_minus,original_v_caret_ea_e_minus);
    [~,maxJ] = Outlier_Detection(length(valid_satellite_numbers),H_G_e,delta_z_min,sigma_rho);
    [~,maxJDot] = Outlier_Detection(length(valid_satellite_numbers),H_G_e,delta_z_dot_min,sigma_rho_dot);
end

% Set initial estimation with all outliers removed
x_est = x_est_m;
% Detect satellites that has outliers for following epochs
for epoch=2:epoch_num
    epoch_pseudo_range_measurements = pseudo_ranges(epoch,:);
    epoch_pseudo_range_rate_measurements = pseudo_range_rates(epoch,:);
    [~,H_G_e,delta_z_min,delta_z_dot_min,r_caret_ea_e_minus,v_caret_ea_e_minus] = GNSS_Least_Squares(times(epoch),epoch_pseudo_range_measurements,epoch_pseudo_range_rate_measurements,satellite_numbers,r_caret_ea_e_minus,v_caret_ea_e_minus);
    [epoch_valid_psuedo_range_indicies,~] = Outlier_Detection(length(satellite_numbers),H_G_e,delta_z_min,sigma_rho);
    [epoch_valid_psuedo_range_rate_indicies,~] = Outlier_Detection(length(satellite_numbers),H_G_e,delta_z_dot_min,sigma_rho_dot);
    valid_pseudo_range_indicies(epoch,:) = epoch_valid_psuedo_range_indicies;
    valid_pseudo_range_rate_indicies(epoch,:) = epoch_valid_psuedo_range_rate_indicies;
end

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