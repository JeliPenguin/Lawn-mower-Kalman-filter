function filtered_heading = new(dr_measurement_data,times)

Define_Constants;

dr_measurement = dr_measurement_data(:,2:end);

% Kalman Filter for Sensor Fusion

% Initialization
tau = tau;  % Time step in seconds
Q_gyro = 0.01;  % Process noise covariance for gyroscope
Q_mag = sigma_mh^2;  % Process noise covariance for magnetometer
R_mag = sigma_mh^2;  % Measurement noise covariance for magnetometer

% Initial state estimate and covariance
x_hat = 0;      % Initial heading estimate
P = 1;          % Initial covariance estimate



% Generate sensor data (for demonstration purposes)
gyro_data = dr_measurement(:,5);  % Replace this with actual gyroscope data
mag_data = dr_measurement(:,6)*deg_to_rad;  % Replace this with actual magnetometer data

% Simulation parameters
num_samples = height(times);

filtered_heading = zeros(num_samples,1);
filtered_heading(1) = mag_data(1);

% Kalman Filter Loop
for k = 2:num_samples
    % Prediction Step (state and covariance prediction)
    x_hat_minus = x_hat + gyro_data(k-1) * tau;
    P_minus = P + Q_gyro;
    
    % Update Step (Kalman gain, measurement residual, and update state and covariance)
    K = P_minus / (P_minus + R_mag);
    residual = mag_data(k) - x_hat_minus;
    x_hat = x_hat_minus + K * residual;
    P = (1 - K) * P_minus;
    
    % Save the filtered heading estimate
    filtered_heading(k) = x_hat*rad_to_deg;

    disp([filtered_heading(k),mag_data(k)*rad_to_deg])
end

% % Plot the results
% figure;
% t = (1:num_samples) * tau;
% subplot(2, 1, 1);
% plot(t, gyro_data, 'b', t, filtered_heading, 'r');
% legend('Gyroscope Data', 'Filtered Heading');
% xlabel('Time (s)');
% ylabel('Angular Rate (deg/s) / Heading (deg)');
% title('Sensor Fusion using Kalman Filter');
% 
% subplot(2, 1, 2);
% plot(t, mag_data, 'b', t, filtered_heading, 'r');
% legend('Magnetometer Data', 'Filtered Heading');
% xlabel('Time (s)');
% ylabel('Magnetic Field Strength / Heading (deg)');
% title('Sensor Fusion using Kalman Filter');

end