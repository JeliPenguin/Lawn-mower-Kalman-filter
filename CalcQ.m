function Q = CalcQ(tau,S_a,S_cf,S_c_phi)
% Calculates system noise covariance matrix Q
% Inputs:
%   tau                 Propagation interval (s)
%   S_a                 Acceleration power spectral density
%   S_cf                Clock frequency PSD
%   S_c_phi             Clock phase PSD
%
% Outputs:
%   Phi                 System noise covariance matrix 8x8

Q = [S_a*tau^3*eye(3)/3,S_a*tau^2*eye(3)/2,zeros(3,1),zeros(3,1)
     S_a*tau^2*eye(3)/2,S_a*tau*eye(3),zeros(3,1),zeros(3,1)
     zeros(1,3),zeros(1,3),S_c_phi*tau+S_cf*tau^3/3,S_cf*tau^2/2
     zeros(1,3),zeros(1,3),S_cf*tau^2/2,S_cf*tau
    ];
end