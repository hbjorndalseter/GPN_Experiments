%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by Andrea Valmorbida and Anese Giovanni
% DII and CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;
addpath(genpath(pwd));


%% Time

[gpsweek, tow, doy, dow] = greg2gps([2024, 1, 31, 0, 0, 0]);


%% Initial settings

% SP3 file (http://navigation-office.esa.int/products/gnss-products/)
filename = 'sp3_files/ESA0MGNFIN_20240310000_01D_05M_ORB.SP3';

% Selected constellation
constID = "GPS";

% Satellite IDs to study
SatID = [1, 15, 30];

% Time vector for interpolation [hours]
interp_time = linspace(0, 3, 1000);


%% Read SP3 file

[sp3, Greg_time] = read_sp3_multiconstellation(filename, constID);


%% Convert SP3 ECEF positions to ECI (stored in sp3.GPS_ECI)

% sp3 = sp3_compute_sc_pos_eci(sp3, Greg_time, constID);


%% Extract per-satellite ECEF positions from SP3

[SAT, SatID] = sp3_get_sc_pos(sp3, SatID);


%% 3-D orbit visualisation with Satellite Scenario

sc           = satelliteScenario;
sc.StartTime = datetime(Greg_time(1, :));
sc.StopTime  = datetime(Greg_time(1, :)) + days(1);

fname = fieldnames(SAT);
% Time axis [s] referenced to first SP3 epoch
t_sp3 = (juliandate(Greg_time) - juliandate(Greg_time(1, :))) * 86400;

for i = 1:length(SatID)
    posTS = timeseries(SAT.(string(fname(i))) * 1e3, t_sp3);  % km → m
    satellite(sc, posTS, "CoordinateFrame", "ecef", "Name", string(fname(i)));
end
satelliteScenarioViewer(sc);


%% Interpolation of SP3 positions (Lagrange, 10 pts)

for i = 1:length(SatID)
    str   = string(fname(i));
    pos_i = SAT.(str);                             % Nx3 [km]
    [orbit_i] = interp_precise_orbits(pos_i, t_sp3, interp_time);
    SAT.(string([char(str), '_interp'])) = orbit_i; % Mx3 [km]
end


%% Plot: SP3 discrete points vs interpolated curve (X component)

interp_time_s = interp_time * 3600;

figure('Name', 'SP3 vs Interpolated Orbits');
for i = 1:length(SatID)
    str      = string(fname(i));
    pos_sp3  = SAT.(str);
    pos_int  = SAT.(string([char(str), '_interp']));

    subplot(length(SatID), 1, i);
    plot(t_sp3,          pos_sp3(:, 1),  'o', 'DisplayName', 'SP3 points'); hold on;
    plot(interp_time_s,  pos_int(:, 1),  '-', 'DisplayName', 'Interpolated');
    xlabel('time [s]');  ylabel('X [km]');
    title(['Satellite ', num2str(SatID(i)), ' — X coordinate']);
    legend;  grid on;
    xlim([0 5000]);

end

%% Save results
save('results_precise_GLONASS.mat', 'SAT', 'SatID', 't_sp3', 'interp_time');


