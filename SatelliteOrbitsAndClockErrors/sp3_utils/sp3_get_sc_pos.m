%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code extracts ECEF positions of desired satellites from the sp3
% structure produced by read_sp3_multiconstellation.
%
% Inputs:   sp3    - structure with sp3.data table:
%                    cols = [Date | GPS_week | GPS_TOW | PRN | x | y | z | clk]
%           SatID  - number(s) of satellite(s) to be investigated
%
% Outputs:  SAT    - structure with ECEF positions [km] per satellite
%                    field name: 'SatID_<N>', value: Nx3 matrix [x y z]
%           SatID  - Available SatID number(s) (unavailable ones removed)
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SAT, SatID] = sp3_get_sc_pos(sp3, SatID)

% Preallocation: indices of unavailable satellites to discard
i_el = [];

for i = 1:length(SatID)

    % ── Check availability of this satellite in the SP3 data ──────────
    % PRN is in column 4 of sp3.data
    if nnz(sp3.data.PRN == SatID(i)) == 0
        disp(['SatID = ', num2str(SatID(i)), ' not found in SP3'])
        i_el = [i_el, i]; %#ok<AGROW>
        continue
    end

    % ── Field name for this satellite ─────────────────────────────────
    str = string(['SatID_', num2str(SatID(i), '%i')]);

    % ── Extract all rows belonging to this PRN, keep x y z (cols 5-7) ─
    mask = sp3.data.PRN == SatID(i);
    GS   = sp3.data(mask, 5:7);          % sub-table with x, y, z [km]

    % ── Store as numeric array in SAT structure ────────────────────────
    SAT.(str) = table2array(GS);          % Nx3 double [km]
end

% ── Remove unavailable satellite IDs from output ──────────────────────
SatID(i_el) = [];

end