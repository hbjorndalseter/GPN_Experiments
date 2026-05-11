%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes derivative of the state vector for the GLONASS
% constellation considering broadcast data
%
% Inputs:   y - state vector
%           acc - vector with pertubing accelerations
%
% Outputs:  dydt - derivative of state vector
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = GLONASS_eq(y,acc)

% Parameters
mu = 3.9860044e11;
J2 = -1.08263e-3;
w_e = 0.7292115e-4;
r_e = 6.378136e6;

...

end