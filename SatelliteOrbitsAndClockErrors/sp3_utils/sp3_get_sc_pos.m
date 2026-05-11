%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code updates the sp3 structure with position of desired
% satellites in ECI reference frame
%
% Inputs:   sp3 - structure with info about sp3 data [GPS_week GPS_TOW PRN x y z]
%           SatID - number(s) of staellite to be investigated
%
% Outputs:  SAT - structure with positions in ECEF for each available 
%                 satellite of SatID
%           SatID - Avilable SatID number(s)
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SAT,SatID] = sp3_get_sc_pos(sp3,SatID)

% preallocation for SatID numbers to discard
i_el = [];


for i = 1 : length(SatID)
    % Control the avilability of SatID data
    if nnz(...)==0 % data not avilable
        disp(['SatID = ', num2str(SatID(i)), ' not found in SP3'])
        i_el = [i_el, i];
        continue
    end
    % Strig with SatID number
    str = string(['SatID_', num2str(SatID(i),'%i')]);
    % Data of SatID [z,y,z]
    GS = ...;
    % Save GS data in SAT structure
    SAT.(str) = table2array(GS);
end

% Discard non available SatID numbers
SatID(i_el) = [];

end