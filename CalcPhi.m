function Phi = CalcPhi(tau)
% Calculates transition matrix Phi
% Inputs:
%   tau                 Propagation interval (s)
%
% Outputs:
%   Phi                 Transition Matrix 8x8

Phi = [eye(3),eye(3)*tau,zeros(3,1),zeros(3,1)
       zeros(3),eye(3),zeros(3,1),zeros(3,1)
       zeros(1,3),zeros(1,3),1,tau
       zeros(1,3),zeros(1,3),0,1];

end