%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare broadcast (RINEX) vs precise (SP3) satellite orbits.
% Computes and plots the 3-D position error of broadcast ephemeris
% with respect to the interpolated SP3 reference.
%
% Coded by Andrea Valmorbida and Anese Giovanni
% DII and CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;
addpath(genpath(pwd));


%% Settings

sp3_file  = 'sp3_files/ESA0MGNFIN_20240310000_01D_05M_ORB.SP3';
brdc_file = 'rinex_obs_files/BRDC00IGS_R_20240310000_01D_MN.rnx';

constID   = "GLONASS";

% Map to RINEX
switch constID
    case "GPS",     rdxID = "GPS";
    case "GALILEO", rdxID = "Galileo";
    case "GLONASS", rdxID = "GLONASS";
    case "BEIDOU",  rdxID = "BeiDou";
end

SatID     = [2, 7, 9];

% Common time vector [hours] for comparison
interp_time = linspace(0, 3, 1000);

if rdxID == "GLONASS"
    interp_time = linspace(0.26, 3.26, 1000); 
end
toe_start   = 259200;                  % [s] GPS TOW at start of day


%% ── PRECISE ORBITS (SP3 + Lagrange interpolation) ────────────────────

[sp3, Greg_time] = read_sp3_multiconstellation(sp3_file, constID);
[SAT_sp3, SatID] = sp3_get_sc_pos(sp3, SatID);

fname_sp3 = fieldnames(SAT_sp3);
t_sp3     = (juliandate(Greg_time) - juliandate(Greg_time(1, :))) * 86400;  % [s]

orbit_precise = cell(length(SatID), 1);
for i = 1:length(SatID)
    str = string(fname_sp3(i));
    orbit_precise{i} = interp_precise_orbits(SAT_sp3.(str), t_sp3, interp_time);
end


%% ── BROADCAST ORBITS (RINEX navigation) ──────────────────────────────

brdc      = rinexread(brdc_file);
info      = rinexinfo(brdc_file);
brdcEphem = brdc.(rdxID);

[SAT_brdc, SatID] = brdc_get_sc_data(brdc, SatID, brdcEphem, rdxID);
fname_brdc = fieldnames(SAT_brdc);

orbit_broadcast = cell(length(SatID), 1);

if rdxID == "GPS" || rdxID == "Galileo" || rdxID == "BeiDou"
    t_vec = toe_start + interp_time * 3600;
    for i = 1:length(SatID)
        str   = string(fname_brdc(i));
        ep    = SAT_brdc.(str).ephTable(1, :);   % first ephemeris epoch
        orb_i = zeros(length(interp_time), 3);
        for j = 1:length(interp_time)
            [x, y, z] = GPS_coordinates(t_vec(j), ...
                ep.Crs,          ep.Delta_n,   ep.M0,      ep.Cuc,  ...
                ep.Eccentricity, ep.Cus,       ep.sqrtA,   ep.Toe,  ...
                ep.Cic,          ep.OMEGA0,    ep.Cis,     ep.i0,   ...
                ep.Crc,          ep.omega,     ep.OMEGA_DOT, ep.IDOT);
            orb_i(j, :) = [x, y, z];
        end
        orbit_broadcast{i} = orb_i;
    end

elseif rdxID == "GLONASS"
    Leap_Seconds = 18;   % GPS - GLONASS time offset [s] for 2024
    for i = 1:length(SatID)
    str    = string(fname_brdc(i));
    ephAll = SAT_brdc.(str).ephTable;

    ep_dt    = ephAll.Properties.RowTimes;
    ep_times = hour(ep_dt)*3600 + minute(ep_dt)*60 + second(ep_dt);

    % Convert all tb times to GPS time
    tb_gps = ep_times + Leap_Seconds;

    % Query times in GPS seconds
    t_query = interp_time * 3600;

    % For each query time, find the nearest tb
    nearest_idx = zeros(length(t_query), 1);
    for j = 1:length(t_query)
        [~, nearest_idx(j)] = min(abs(tb_gps - t_query(j)));
    end

    % Identify segment boundaries (where the chosen ephemeris row changes)
    seg_starts = [1; find(diff(nearest_idx) ~= 0) + 1];
    seg_ends   = [seg_starts(2:end) - 1; length(t_query)];

    x_all = zeros(length(t_query), 1);
    y_all = zeros(length(t_query), 1);
    z_all = zeros(length(t_query), 1);

    for s = 1:length(seg_starts)
        idx_seg = seg_starts(s):seg_ends(s);
        eph_idx = nearest_idx(idx_seg(1));
        ep      = ephAll(eph_idx, :);
        tb_s    = tb_gps(eph_idx);

        t_rel = t_query(idx_seg) - tb_s;

        [x, y, z] = GLONASS_coordinates(ep, t_rel);
        x_all(idx_seg) = x;
        y_all(idx_seg) = y;
        z_all(idx_seg) = z;
    end

    fprintf('SV%d: %d ephemeris segments used over %.1f–%.1f hr\n', ...
            SatID(i), length(seg_starts), interp_time(1), interp_time(end));

    orbit_broadcast{i} = [x_all, y_all, z_all];
    end
end


%% ── POSITION ERROR: broadcast − precise ──────────────────────────────

t_hr = interp_time;   % [hours] common time axis

figure('Name', 'Broadcast vs Precise Orbit Error');
for i = 1:length(SatID)
    err_3d = orbit_broadcast{i} - orbit_precise{i};   % [km]
    err_norm = vecnorm(err_3d, 2, 2) * 1e3;           % [m]

    subplot(length(SatID), 1, i);
    plot(t_hr, err_norm, 'LineWidth', 1.5);
    xlabel('time [hr]');
    ylabel('3-D error [m]');
    title(['Satellite ', num2str(SatID(i)), ...
           ' — Broadcast vs Precise orbit error']);
    grid on;
end

figure('Name', 'Component Errors');
comp_labels = {'X', 'Y', 'Z'};
for i = 1:length(SatID)
    err_3d = (orbit_broadcast{i} - orbit_precise{i}) * 1e3;   % [m]
    for c = 1:3
        subplot(length(SatID), 3, (i-1)*3 + c);
        plot(t_hr, err_3d(:, c), 'LineWidth', 1.2);
        xlabel('time [hr]');
        ylabel([comp_labels{c}, ' error [m]']);
        title(['SV', num2str(SatID(i)), ' Δ', comp_labels{c}]);
        grid on;
    end
end

figure('Name', 'Single segment zoom');
err_3d = vecnorm((orbit_broadcast{1} - orbit_precise{1}) * 1e3, 2, 2);
plot(interp_time, err_3d, 'LineWidth', 1.2);
xlim([0.5, 1.0]);
xlabel('time [hr]'); ylabel('3-D error [m]');
title('SV2 — Single segment zoom');
grid on;

figure('Name', 'Broadcast vs Precise — zoom');
subplot(2,1,1);
plot(interp_time, orbit_broadcast{1}(:,1), 'LineWidth', 1.2);
xlim([0.5, 1.0]);
xlabel('time [hr]'); ylabel('X [km]');
title('SV2 — Broadcast X (zoom)');
grid on;

subplot(2,1,2);
plot(interp_time, orbit_precise{1}(:,1), 'LineWidth', 1.2);
xlim([0.5, 1.0]);
xlabel('time [hr]'); ylabel('X [km]');
title('SV2 — Precise X (zoom)');
grid on;

% ── Summary statistics ─────────────────────────────────────────────────
results_table = table();
results_table.SatID = SatID(:);
results_table.RMS_m  = zeros(length(SatID), 1);
results_table.Max_m  = zeros(length(SatID), 1);
results_table.Mean_m = zeros(length(SatID), 1);

fprintf('\n%-6s  %10s  %10s  %10s\n', 'SatID', 'RMS [m]', 'Max [m]', 'Mean [m]');
for i = 1:length(SatID)
    err = vecnorm((orbit_broadcast{i} - orbit_precise{i}) * 1e3, 2, 2);
    results_table.RMS_m(i)  = rms(err);
    results_table.Max_m(i)  = max(err);
    results_table.Mean_m(i) = mean(err);
    fprintf('%-6d  %10.3f  %10.3f  %10.3f\n', ...
            SatID(i), results_table.RMS_m(i), results_table.Max_m(i), results_table.Mean_m(i));
end

% ── Save results ────────────────────────────────────────────────────────
if ~exist('results', 'dir'), mkdir('results'); end

% Save workspace variables
save(fullfile('results', ['compare_', char(constID), '.mat']), ...
     'orbit_precise', 'orbit_broadcast', 'SatID', 'interp_time', 'results_table');

% Save summary table as CSV
writetable(results_table, fullfile('results', ['compare_', char(constID), '.csv']));

% Save summary as readable text file
fid = fopen(fullfile('results', ['compare_', char(constID), '.txt']), 'w');
fprintf(fid, 'Constellation: %s\n', constID);
fprintf(fid, 'Time window:   %.1f hours\n\n', interp_time(end));
fprintf(fid, '%-6s  %10s  %10s  %10s\n', 'SatID', 'RMS [m]', 'Max [m]', 'Mean [m]');
for i = 1:length(SatID)
    fprintf(fid, '%-6d  %10.3f  %10.3f  %10.3f\n', ...
            SatID(i), results_table.RMS_m(i), results_table.Max_m(i), results_table.Mean_m(i));
end
fclose(fid);

fprintf('\nResults saved to /results/compare_%s (.mat, .csv, .txt)\n', constID);