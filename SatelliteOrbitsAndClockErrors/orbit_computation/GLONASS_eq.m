%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes derivative of the state vector for the GLONASS
% constellation considering broadcast data
%
% Inputs:   y   - state vector [x; vx; y; vy; z; vz]  (SI units: m, m/s)
%           acc - vector with perturbing accelerations [ax; ay; az] (m/s²)
%
% Outputs:  dydt - derivative of state vector
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = GLONASS_eq(y, acc)
% PZ-90 constants in SI (meters, seconds)
mu  = 3.9860044e14;    % [m³/s²]
J2  = 1.08263e-3;      % [-]
w_e = 7.2921151467e-5; % [rad/s]
r_e = 6378136.0;       % [m]

r = sqrt(y(1)^2 + y(3)^2 + y(5)^2);

dydt = zeros(6,1);

dydt(1) = y(2);
dydt(2) = - mu*y(1)/r^3 ...
          + (3/2)*J2*(mu*r_e^2/r^5)*y(1)*(1 - 5*y(5)^2/r^2) ...
          + w_e^2*y(1) + 2*w_e*y(4) ...
          + acc(1);

dydt(3) = y(4);
dydt(4) = - mu*y(3)/r^3 ...
          + (3/2)*J2*(mu*r_e^2/r^5)*y(3)*(1 - 5*y(5)^2/r^2) ...
          + w_e^2*y(3) - 2*w_e*y(2) ...
          + acc(2);

dydt(5) = y(6);
dydt(6) = - mu*y(5)/r^3 ...
          + (3/2)*J2*(mu*r_e^2/r^5)*y(5)*(3 - 5*y(5)^2/r^2) ...
          + acc(3);
end