%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Builds a SAT structure with broadcast ephemeris data for the
% requested satellites of a given constellation.
%
% Inputs:
%   brdc      - full structure returned by rinexread (not used directly
%               here, kept for interface compatibility)
%   SatID     - vector of requested satellite PRN numbers
%   brdcEphem - timetable/table for the selected constellation,
%               i.e. brdc.(constID) from rinexread output
%   constID   - constellation string: "GPS", "GLONASS", "Galileo", "BeiDou"
%
% Outputs:
%   SAT   - structure with one field per available satellite.
%           Field name : 'SatID_<N>'
%           Field value: struct with .ID and .ephTable (rows of brdcEphem)
%   SatID - updated vector containing only the available satellite IDs
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SAT, SatID] = brdc_get_sc_data(brdc, SatID, brdcEphem, constID) %#ok<INUSL>
% brdcEphem is the constellation timetable (e.g. brdc.GPS)
ephTable = brdcEphem;
% ── Find which requested satellites are actually present ──────────────
availableSatIDs = unique(ephTable.SatelliteID);
validSatID      = intersect(SatID, availableSatIDs);
if isempty(validSatID)
    warning('brdc_get_sc_data: None of the requested SatIDs found in %s ephemeris.', constID);
    SAT   = struct();
    SatID = [];
    return;
end
SatID = validSatID;    % keep only satellites that exist in the file
% ── Build SAT structure ───────────────────────────────────────────────
SAT = struct();
for i = 1:length(SatID)
    mask = ephTable.SatelliteID == SatID(i);
    str  = string(['SatID_', num2str(SatID(i), '%i')]);
    SAT.(str).ID       = SatID(i);
    SAT.(str).ephTable = ephTable(mask, :);   % all rows for this SV
end
end