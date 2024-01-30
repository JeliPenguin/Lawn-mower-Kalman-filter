function avg_speed = Calc_Avg_Speed(wheel_speed_measurements)
% Calculates average speed of vehicle based on speed measurements for each wheel
% Inputs:
%   wheel_speed_measurements        Array of speed measurements for each wheel
% Output:
%   avg_speed                       Calculated average speed
avg_speed = mean(wheel_speed_measurements);
end