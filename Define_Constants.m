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

sigma_rho = sqrt(sigma_rho_assumed^2+sigma_space^2+sigma_res_ion^2+sigma_res_trop^2+sigma_code_multipath^2); % root sum of squares of all relevant error sources
sigma_rho_dot = sqrt(sigma_rho_dot_assumed^2+sigma_rho_dot_multipath^2); % root sum of squares of all relevant error sources

sigma_p = 5; % Measurement error standard deviation (m)
S_c_phi = 0.01; % Clock phase PSD
S_cf = 0.04; % Clock frequency PSD

T = 3; % Outlier detection threshold

tau = 0.5; % Propagation interval
S_a = 0.01; % Acceleration power spectral density

sigma_Gr = 5; % GNSS position measurements error standard deviation
sigma_Gv = 0.02; % GNSS velocity measurements have an error standard deviation
S_DR = 0.2; % DR velocity error power spectral density

% Ends