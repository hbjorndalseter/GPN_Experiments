%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code interpolates satellite positions from an SP3 file data
%
% Inputs:   pos - matrix with satellite coordinates from SP3 file [km]
%           time - time span vector [s]
%           interp_time - vector of time where compute the interpolation [hours]
%
% Outputs:  orbit - interpolated coordinates of satellite orbit [km]
%            
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [orbit] = interp_precise_orbits(pos,time,interp_time)

...
