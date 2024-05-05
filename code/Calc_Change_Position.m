function [delta_r_rear] = Calc_Change_Position(wheel_speed_measurements, gyro_measurement, compass_measurement)
% Calculates velocity of vehicle based on speed measurements for each wheel
% Inputs:
%   wheel_speed_measurements        Array of speed measurements for each wheel (m/s)
%   gyro_measurement                 gyroscope angular rate measurements (rad/s)
%   compass_measurement              magnetic compass measurements (degrees)
% Output:
%   velocity                       Calculated velocity

Define_Constants;

% Define inputs
psi_nb = compass_measurement * deg_to_rad;
psi_dot_nb = gyro_measurement;

% Calculate driving (rear) wheel-frame average speed
v_erL = wheel_speed_measurements(1);
v_erR = wheel_speed_measurements(2);
v_er = 0.5 * (v_erL + v_erR);

% Calculate change in position of rear wheels
delta_r_rear_temp = [cos(psi_nb) - 0.5 * psi_dot_nb * tau * sin(psi_nb); ...
    sin(psi_nb) + 0.5 * psi_dot_nb * tau * cos(psi_nb)] * v_er * tau; 

delta_r_rear = delta_r_rear_temp + [cos(psi_nb), -sin(psi_nb); sin(psi_nb), cos(psi_nb)] * [l_bry; -l_brx] * psi_dot_nb * tau;

end