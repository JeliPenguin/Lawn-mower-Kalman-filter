function [valid_indicies, maxJ] = Outlier_Detection(number_of_satellites,H_G_e,delta_z_min,sigma)
% Detects outliers within given data
%
% Inputs:
%   number_of_satellites         Total count of number of satellites giving 
%                                the measurements
%   H_G_e                        Measurement Matrix           
%   delta_z_min                  Innovation vector  
%   sigma                        Measurement error standard deviation         
%
% Outputs:
%   valid_indicies               Array of (1,number_of_satellites)
%                                indicating whether the measurement of the
%                                satellite is an outlier (0) or inliner (1)
%   maxJ                         Index of satellite giving the highest
%                                normalized residual value larger than
%                                threshold

Define_Constants;

% Compute residuals vector
v = (H_G_e * inv(H_G_e.' * H_G_e) * H_G_e.' - eye(number_of_satellites)) * delta_z_min; 

% Compute residuals covariance matrix
C_v = (eye(number_of_satellites) - H_G_e * inv(H_G_e.' * H_G_e) * H_G_e.') * sigma^2;

% Setup
valid_indicies = zeros(1,number_of_satellites);
maxResidual = -inf;
maxJ = 0;

% Check for outliers
for j=1:number_of_satellites
    normalized_residual = norm(v(j))/sqrt(C_v(j,j));
    if normalized_residual > T
        if normalized_residual > maxResidual
            maxJ = j;
            maxResidual = normalized_residual;
        end
    else
        valid_indicies(j) = 1;
    end
end
end

