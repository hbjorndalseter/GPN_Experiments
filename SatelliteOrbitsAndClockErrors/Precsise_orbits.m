%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by Andrea Valmorbida and Anese Giovanni
% DII and CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% workspace setting

close all; clear; clc;

addpath(genpath('.\'));


%% Time

[gpsweek, tow, doy, dow] = greg2gps([2024,1,31,0,0,0]);


%% Initial settings

% SP3 files available at http://navigation-office.esa.int/products/gnss-products/
filename = 'ESA0MGNFIN_20240310000_01D_05M_ORB.SP3';

% Selected constellation
constID = "GPS";

% Satellite ID (choose satellite ID to study)
SatID = [1,15,30];

% Time to interpolate
interp_time = linspace(0,3,1000);  % [hours]


%% Read data from SP3 file

% Read SP3 file
[sp3, Greg_time] = ...;


%% Data from SP3 file in ECI reference frame

sp3 = sp3_compute_sc_pos_eci(...);


%% get satellite coordinates from sp3 data format

[SAT,SatID] = sp3_get_sc_pos(...);


%% Visualize orbits
sc = satelliteScenario;
sc.StartTime = datetime(Greg_time(1,:));
sc.StopTime = datetime(Greg_time(1,:))+days(1);
% sat
fname = fieldnames(SAT);
t = (juliandate(Greg_time)-juliandate(Greg_time(1,:)))*86400; % time
for i = 1 : length(SatID)
    postimeseries = timeseries(SAT.(string(fname(i)))*1e3, t);
    satellite(sc,postimeseries,"CoordinateFrame","ecef","Name",fname(i));
end
satelliteScenarioViewer(sc);


%% Interpolation of SP3

% create function that interpolates sat pos based on time to interp, GS_ID
for i = 1 : length(SatID)
    [orbit] = interp_precise_orbits(...);
    ...
end


%% Compare SP3 data and interpolated data 

figure();

...
