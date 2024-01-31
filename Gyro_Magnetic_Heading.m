function gyro_mag_heading_solution = Gyro_Magnetic_Heading(dr_measurement_data,times)
% Calculates gyro smoothed magnetic heading
% Inputs:
%   dr_measurement_data                 Array of raw dead reckoning measurement data
%   times                               Array of time value for each epoch
%
% Outputs:
%   gyro_mag_heading_solution           Array of gyro smoothed magnetic heading solutions

Define_Constants;

dr_measurement = dr_measurement_data(:,2:end);

% Total count of epochs
epoch_num = height(times); 

% Initial heading solution
smoothedHeading = dr_measurement(1,6)*deg_to_rad;

gyro_mag_heading_solution = zeros(epoch_num,1);
gyro_mag_heading_solution(1) = dr_measurement(1,6);

W_m = sigma_gyro*tau/sigma_mh;

for i = 2:epoch_num
    gyro_measurement = dr_measurement(i,5); % rad/s
    magnetic_heading_measurement = dr_measurement(i,6)*deg_to_rad; %rad
    
    smoothedHeading = W_m*magnetic_heading_measurement+(1-W_m)*(smoothedHeading+gyro_measurement*tau);
    gyro_mag_heading_solution(i) = smoothedHeading*rad_to_deg;
end


% Writing solution to csv file
outputTable = array2table([dr_measurement(:,6),gyro_mag_heading_solution]);
writetable(outputTable,"Solutions/Gyro_Smoothed_Mag_Heading_Solution.csv",'WriteVariableNames',0)

end