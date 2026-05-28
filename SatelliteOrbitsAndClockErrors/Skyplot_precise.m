%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by Andrea Valmorbida and Anese Giovanni
% DII and CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;
addpath(genpath(pwd));


%% Settings

% SP3 file (http://navigation-office.esa.int/products/gnss-products/)
filename = 'sp3_files/ESA0MGNFIN_20240310000_01D_05M_ORB.SP3';

% Selected constellation
constID = "GLONASS";

% Satellite IDs to include in the skyplot
% SatID = [1, 2, 3, 4, 5, 6, 14, 15, 16, 28, 29, 30];
SatID = 1:32;


%% Read SP3 file

[sp3, Greg_time] = read_sp3_multiconstellation(filename, constID);


%% Extract per-satellite ECEF positions

[SAT, SatID] = sp3_get_sc_pos(sp3, SatID);


%% Azimuth and elevation computation

% Masking angle [deg]
maskAngle = 15;

% Time axis [s] from first SP3 epoch
t = (juliandate(Greg_time) - juliandate(Greg_time(1, :))) * 86400;

% Receiver position (University of Padova, Italy): [lat°, lon°, alt m]
RecPos_lla_ref = [45.4108534036746, 11.8917094896601, 70.0413117404738];

% Pre-allocate azimuth / elevation / visibility matrices
%   rows = time epochs,  cols = satellites
Az_matrix  = NaN(numel(t), numel(SatID));
El_matrix  = Az_matrix;
Vis_matrix = Az_matrix;

fname = fieldnames(SAT);

% Convert receiver LLA → ECEF [km]
rec_ecef_km = lla2ecef_km(RecPos_lla_ref);

for i = 1:length(SatID)
    sat_ecef_km = SAT.(string(fname(i)));   % Nt×3 [km]

    [az, el, vis] = compute_azel(sat_ecef_km, rec_ecef_km, maskAngle);

    Az_matrix(:, i)  = az;
    El_matrix(:, i)  = el;
    Vis_matrix(:, i) = vis;
end

% Replace sub-zero elevation with NaN (satellite below horizon)
El_matrix(El_matrix < 0) = NaN;


%% Animated SkyPlot

figure('Name', 'SkyPlot');
sp = skyplot([], [], 'MaskElevation', maskAngle);
set(sp, 'LabelFontSize', 20);

satLabels = arrayfun(@(id) ['G', num2str(id, '%02d')], SatID, 'UniformOutput', false);

for idx = 1:size(Az_matrix, 1)

    az_now = Az_matrix(idx, :);
    el_now = El_matrix(idx, :);

    % Only plot satellites that are visible (non-NaN elevation)
    vis_mask = ~isnan(el_now);

    if any(vis_mask)
        sp.AzimuthData   = az_now(vis_mask);
        sp.ElevationData = el_now(vis_mask);
        sp.LabelData     = satLabels(vis_mask);
    end

    drawnow limitrate
    pause(0.05);
end


%% Number of visible satellites over time

nSatView = sum(Vis_matrix, 2);   % sum across columns (satellites)

figure('Name', 'Visible Satellites');
plot(t / 3600, nSatView, 'LineWidth', 1.5);
xlabel('time [hr]');
ylabel('n. of visible s/c');
title('Number of GPS satellites above mask angle');
grid on;


% =========================================================================
% ── Local helper functions ────────────────────────────────────────────────
% =========================================================================

function rec_ecef = lla2ecef_km(lla)
% Convert receiver position from geodetic LLA [deg, deg, m] to ECEF [km].
% Uses WGS-84 ellipsoid parameters.

a   = 6378.137;          % [km]  semi-major axis
f   = 1 / 298.257223563; % [-]   flattening
e2  = 2*f - f^2;         % eccentricity squared

lat = deg2rad(lla(1));
lon = deg2rad(lla(2));
h   = lla(3) / 1e3;      % m → km

N = a / sqrt(1 - e2 * sin(lat)^2);   % prime vertical radius of curvature [km]

rec_ecef = [(N + h)           * cos(lat) * cos(lon); ...
            (N + h)           * cos(lat) * sin(lon); ...
            (N * (1 - e2) + h) * sin(lat)];
end


function [az, el, vis] = compute_azel(sat_ecef_km, rec_ecef_km, maskAngle)
% Compute azimuth [deg], elevation [deg], and visibility flag for a
% satellite whose ECEF positions are given over time.
%
%   sat_ecef_km  - Nt×3 satellite ECEF positions [km]
%   rec_ecef_km  - 3×1 receiver ECEF position    [km]
%   maskAngle    - elevation mask angle           [deg]

Nt = size(sat_ecef_km, 1);
az  = zeros(Nt, 1);
el  = zeros(Nt, 1);
vis = false(Nt, 1);

% Receiver geodetic coordinates (needed for ENU rotation)
rec = rec_ecef_km(:);           % 3×1
[lat, lon] = ecef2latlon(rec);  % [rad]

sinLat = sin(lat);  cosLat = cos(lat);
sinLon = sin(lon);  cosLon = cos(lon);

% Rotation matrix ECEF → ENU
R_enu = [-sinLon,              cosLon,             0;
         -sinLat*cosLon,  -sinLat*sinLon,  cosLat;
          cosLat*cosLon,   cosLat*sinLon,  sinLat];

for k = 1:Nt
    dR = sat_ecef_km(k, :)' - rec;   % line-of-sight vector [km]

    enu = R_enu * dR;                 % East, North, Up components
    E   = enu(1);
    N   = enu(2);
    U   = enu(3);

    el(k) = rad2deg(atan2(U, sqrt(E^2 + N^2)));
    az(k) = mod(rad2deg(atan2(E, N)), 360);   % [0, 360)

    vis(k) = el(k) >= maskAngle;
end
end


function [lat, lon] = ecef2latlon(ecef)
% Approximate conversion ECEF [km] → geodetic latitude and longitude [rad].
x = ecef(1);  y = ecef(2);  z = ecef(3);
lon = atan2(y, x);
p   = sqrt(x^2 + y^2);
lat = atan2(z, p * (1 - 0.00669437999014));   % first approximation (spheroid)
% Bowring iteration for better accuracy
a  = 6378.137;
e2 = 0.00669437999014;
for iter = 1:5
    N   = a / sqrt(1 - e2 * sin(lat)^2);
    lat = atan2(z + e2 * N * sin(lat), p);
end
end