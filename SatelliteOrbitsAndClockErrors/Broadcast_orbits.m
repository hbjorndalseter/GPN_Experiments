%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by Andrea Valmorbida and Anese Giovanni
% DII and CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;
addpath(genpath(pwd));


%% Initial settings

% RINEX navigation file (https://cddis.nasa.gov/archive/gnss/data/daily/)
filename = 'rinex_obs_files/BRDC00IGS_R_20240310000_01D_MN.rnx';

% Selected constellation  ("GPS" | "GLONASS" | "Galileo" | "BeiDou")
constID = "GLONASS";

% Map to rinexread field name (different capitalization)
switch constID
    case "GPS",     rdxID = "GPS";
    case "GALILEO", rdxID = "Galileo";
    case "GLONASS", rdxID = "GLONASS";
    case "BEIDOU",  rdxID = "BeiDou";
end

% Satellite IDs to study
SatID = [1, 10, 15];

% Time vector [hours] and GPS TOW reference epoch
vector_time = linspace(0, 3, 1000);   % [hours]
toe_start   = 259200;                  % [s] start of day in GPS TOW


%%  Read broadcast RINEX navigation file

brdc = rinexread(filename);            % returns struct with .GPS, .GLONASS, etc.
info = rinexinfo(filename);            % header metadata (leap seconds, etc.)

% Ephemeris table for the chosen constellation
brdcEphem = brdc.(rdxID);

%% Extract data for the selected satellites

[SAT, SatID] = brdc_get_sc_data(brdc, SatID, brdcEphem, rdxID);


%% Orbit computation

orbit = zeros(length(vector_time), 3);
fname = fieldnames(SAT);

% ── GPS (and Galileo / BeiDou: same Keplerian algorithm) ──────────────
if rdxID == "GPS" || rdxID == "Galileo" || rdxID == "BeiDou"

    for i = 1:length(SatID)
        str   = string(fname(i));
        ephem = SAT.(str).ephTable;

        % Choose the ephemeris epoch closest to the middle of the interval
        % (in practice the first available row is used here)
        ep = ephem(1, :);

        % GPS time vector [s]
        t_vec = toe_start + vector_time * 3600;

        for j = 1:length(vector_time)
            [x, y, z] = GPS_coordinates(t_vec(j), ...
                ep.Crs,          ep.Delta_n,   ep.M0,      ep.Cuc,  ...
                ep.Eccentricity, ep.Cus,       ep.sqrtA,   ep.Toe,  ...
                ep.Cic,          ep.OMEGA0,    ep.Cis,     ep.i0,   ...
                ep.Crc,          ep.omega,     ep.OMEGA_DOT, ep.IDOT);
            orbit(j, :) = [x, y, z];   % [km]
        end

        % Store broadcast orbit in SAT structure
        SAT.(string([char(str), '_brdc'])) = orbit;
    end

% ── GLONASS (numerical integration of equations of motion) ────────────
elseif rdxID == "GLONASS"
    for i = 1:length(SatID)
        str    = string(fname_brdc(i));
        ephAll = SAT_brdc.(str).ephTable;

        % Get all epoch times in seconds from UTC midnight
        ep_dt    = ephAll.Properties.RowTimes;
        ep_times = hour(ep_dt)*3600 + minute(ep_dt)*60 + second(ep_dt);

        % Select the epoch closest to the start of the interpolation window
        [~, idx] = min(abs(ep_times - interp_time(1)*3600));
        ep   = ephAll(idx, :);
        tb_s = ep_times(idx);

        % Time relative to ephemeris epoch [s] — t=0 means "at tb"
        t_rel = interp_time * 3600 - tb_s;

        [x, y, z] = GLONASS_coordinates(ep, t_rel);
        orbit_broadcast{i} = [x, y, z];
    end
end


%% Visualise broadcast orbits (X component vs time)

figure('Name', ['Broadcast Orbits — ', char(constID)]);
t_plot = vector_time * 3600;   % [s]

for i = 1:length(SatID)
    str       = string(fname(i));
    orbit_plt = SAT.(string([char(str), '_brdc']));

    subplot(length(SatID), 1, i);
    plot(t_plot, orbit_plt(:, 1), '-');
    xlabel('time [s]');  ylabel('X [km]');
    title(['Satellite ', num2str(SatID(i)), ' — Broadcast orbit X']);
    grid on;
    % xlim([0, 5000]);
end

%% Save results
save('results_broadcast_GPS.mat', 'SAT', 'SatID', 'vector_time');