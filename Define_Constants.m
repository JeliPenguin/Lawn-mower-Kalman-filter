%Define_Constants
%SCRIPT This defines a number of constants for your use
% Created 16/11/16 by Paul Groves

% Copyright 2016, Paul Groves
% License: BSD; see license.txt for details

% Constants
deg_to_rad = 0.01745329252; % Degrees to radians conversion factor
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor
c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]);
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity

sigma_r = 10; % Initial position standard deviation, assuming same for all directions
sigma_v = 0.1; % Initial velocity standard deviation, assuming same for all directions
sigma_co = 100000; % Initial clock offset standard deviation
sigma_cd = 200; % Reciever clock drift standard deviation

sigma_space = 1; % Signal in space error std
sigma_res_ion = 2; % Residual ionosphere error standard deviation at zenith
sigma_res_trop = 0.2; % Residual troposphere error standard deviation at zenith
sigma_code_multipath = 2; % Code tracking and multipath error std
sigma_rho_dot_multipath = 0.02; % Range rate tracking and multipath error std

sigma_rho_assumed = 10; % Assumed Psuedo-range measurement error std
sigma_rho_dot_assumed = 0.05; % Assumed Pseudo-range rate measurement error std

sigma_rho = sigma_rho_assumed;
sigma_rho_dot = sigma_rho_dot_assumed;

%% TO DISCUSS - suggest we use the above sigma_rho not these, based on the results
% sigma_rho = sqrt(sigma_space^2+sigma_res_ion^2+sigma_res_trop^2+sigma_code_multipath^2); % root sum of squares of all relevant error sources
% sigma_rho_dot = sqrt(sigma_rho_dot_multipath^2); % root sum of squares of all relevant error sources

S_c_phi = 0.01; % Clock phase PSD
S_cf = 0.04; % Clock frequency PSD

T = 6; % Outlier detection threshold

tau = 0.5; % Propagation interval
S_a = 0.01; % Acceleration power spectral density

sigma_Gr = 5; % GNSS position measurements error standard deviation
sigma_Gv = 0.02; % GNSS velocity measurements have an error standard deviation
S_DR = 0.2; % DR velocity error power spectral density


micro_g_to_meters_per_second_squared = 9.80665E-6;

% Loosely coupled INS/GNSS Kalman filter parameters
S_rg = (0.02 * deg_to_rad / 60)^2; % Gyro noise PSD (deg^2 per hour, converted to rad^2/s)
% Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)                
LC_KF_config.accel_noise_PSD = (200 *...
    micro_g_to_meters_per_second_squared)^2;
% Accelerometer bias random walk PSD (m^2 s^-5)
LC_KF_config.accel_bias_PSD = 1.0E-7;

S_bgd = 2.0E-12; % Gyro bias random walk PSD (rad^2 s^-3)
% Position measurement noise SD per axis (m)
LC_KF_config.pos_meas_SD = 2.5;
% Velocity measurement noise SD per axis (m/s)
LC_KF_config.vel_meas_SD = 0.05;

sigma_mh = 4*deg_to_rad; % Magnetic heading error standard deviation
sigma_gyro = 0.0001; % Gyroscope angular rate error standard deviation

% Ends