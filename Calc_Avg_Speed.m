function avg_speed = Calc_Avg_Speed(wheel_speed_measurements,gyro_measurement)
% Calculates average speed of vehicle based on speed measurements for each wheel
% Inputs:
%   wheel_speed_measurements        Array of speed measurements for each wheel (m/s)
%   gyro_measurment                 gyroscope angular rate measurements (rad/s)
% Output:
%   avg_speed                       Calculated average speed

%% TODO Apply compensation based on turn rate from gyro measurements, or averaging speeds of left and right wheels
avg_speed = mean(wheel_speed_measurements);
end