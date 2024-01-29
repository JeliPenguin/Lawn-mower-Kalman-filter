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

sigma_r = [10;10;10]; % Initial position standard deviation
sigma_v = [0.1;0.1;0.1]; % Initial velocity standard deviation
sigma_co = 10; % Initial clock offset standard deviation
sigma_cd = 0.1; % Initial clock drift standard deviation
sigma_rho = 10; % Psuedo-range measurement error std
sigma_rho_dot = 0.05; % Pseudo-range rate measurement error std
sigma_p = 5; % Measurement error standard deviation (m)
S_c_phi = 0.01; % Clock phase PSD
S_cf = 0.04; % Clock frequency PSD
T = 6; % Outlier detection threshold
tau = 1; % Propagation interval
S_a = 5; % Acceleration power spectral density

% Ends