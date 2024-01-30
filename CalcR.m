function R = CalcR(sigma_rho,sigma_rho_dot)
% Calculate measurement noise covariance matrix R
% Inputs:
%   sigma_rho                    Psuedo-range measurement error std
%   sigma_rho_dot                Psuedo-range rate measurement error std
%
% Outputs:
%   R                 Measurement noise covariance matrix 8x8

R = [eye(8,8) * sigma_rho^2,zeros(8,8)
     zeros(8,8),eye(8,8)*sigma_rho_dot^2];
end