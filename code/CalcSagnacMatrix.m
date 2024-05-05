function C_I_e = CalcSagnacMatrix(range)
% Calculates Sagnac matrix given the range prediction
% Inputs:
%   range                   Single range prediction value
%
% Outputs:
%   C_I_e                   Calculated Sagnac matrix
    Define_Constants;
    rep = omega_ie * range/c;
    C_I_e = [1,rep,0;-rep,1,0;0,0,1];
end