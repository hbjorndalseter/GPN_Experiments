%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the positions in ECEF reference frame for GLONASS
% constellation starting from broadcast data
%
% Inputs:   ephem - structure with broadcast ephemeris data
%           interp_time - time vector where user wants to compute positions
%
% Outputs:  x - x coordinates in ECEF [km]
%           y - y coordinates in ECEF [km]
%           z - z coordinates in ECEF [km]
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,z] = GLONASS_coordinates(ephem, interp_time)

% Initial data
pos_0 = [ephem.PositionX,ephem.PositionY,ephem.PositionZ]';
vel_0 = [ephem.VelocityX,ephem.VelocityY,ephem.VelocityZ]';
y0 = [pos_0(1), vel_0(1), pos_0(2), vel_0(2), pos_0(3), vel_0(3)]';

% perturbation vector
acc_p = [ephem.AccelerationX,ephem.AccelerationY,ephem.AccelerationZ]';

% Integration of differential equations
[t,y_out] = ...;

% Positions
x = y_out(:,1);
y = y_out(:,3);
z = y_out(:,5);

