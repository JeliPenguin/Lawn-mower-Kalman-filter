function gyro_mag_heading_solution = Gyro_Magnetometer_KF(dr_measurement_data,times)
% Calculates an gyro-smoothed magnetic heading solution using Kalman
% filtering
% Inputs:
%   dr_measurement_data                      Array of raw dead reckoning measurement data
%   times                                    Array of time value for each epoch
%
% Outputs:
%   gyro_mag_heading_solution                Calculated Heading Solution
%                                               Col 1: Heading (deg)

Define_Constants;

dr_measurement = dr_measurement_data(:,2:end);

% Initialization
Q = 0.01;  % Process noise covariance for gyroscope
R_k = sigma_mh^2;  % Measurement noise covariance for magnetometer

% Initial state estimate and covariance
x_hat = 0;      % Initial heading estimate
P = 1;          % Initial covariance estimate

gyro_data = dr_measurement(:,5); 
mag_data = dr_measurement(:,6)*deg_to_rad;

% Simulation parameters
epoch_num = height(times);

gyro_mag_heading_solution = zeros(epoch_num,1);
gyro_mag_heading_solution(1) = mag_data(1);

% Kalman Filter Loop
for k = 2:epoch_num
    % Prediction Step (state and covariance prediction)
    x_hat_minus = x_hat + gyro_data(k-1) * tau;
    P_minus = P + Q;
    
    % Update Step (Kalman gain, measurement residual, and update state and covariance)
    K = P_minus * inv(P_minus + R_k);
    delta_z = mag_data(k) - x_hat_minus;
    x_hat = x_hat_minus + K * delta_z;
    P = (1 - K) * P_minus;
    
    % Save the filtered heading estimate
    gyro_mag_heading_solution(k) = x_hat*rad_to_deg;
end

t = (1:epoch_num) * tau;
f = figure();
plot(t, mag_data*rad_to_deg, 'b', t, gyro_mag_heading_solution, 'r');
legend('Magnetometer Data', 'Filtered Heading');
xlabel('Time (s)');
ylabel('Heading (deg)');
title('Gyro Smoothed magnetic heading');
saveas(f,"Figures/Heading/smoothedComparison.png")

% Writing solution to csv file
outputTable = array2table(gyro_mag_heading_solution);
writetable(outputTable,"Solutions/gyro_smoothed_heading_solution.csv",'WriteVariableNames',0)

end