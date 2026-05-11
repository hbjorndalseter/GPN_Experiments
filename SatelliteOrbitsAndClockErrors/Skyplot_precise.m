%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by Andrea Valmorbida and Anese Giovanni
% DII and CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% workspace setting
close all; clear; clc;

addpath(genpath('.\'));


%%
% file available at http://navigation-office.esa.int/products/gnss-products/

filename = 'ESA0MGNFIN_20240310000_01D_05M_ORB.SP3';

% Selected constellation
constID = "GPS";

% Satellite ID (choose satellite ID to study)
SatID = [1,2,3,4,5,6,14,15,16,28,29,30];
% SatID = 1:1:32;


%% Read data

% Read SP3 file
[sp3, Greg_time] = ...;


%% get satellite coordinates from sp3 data format

[SAT,SatID] = ...;


%% compute azimuth and elevatoin

% set Masking Angle
maskAngle = 15;     % [deg]
t = (juliandate(Greg_time)-juliandate(Greg_time(1,:)))*86400; % s

RecPos_lla_ref = [45.4108534036746, 11.8917094896601, 70.0413117404738];

Az_matrix = NaN(numel(t),numel(SatID));
El_matrix = Az_matrix;
Vis_matrix = Az_matrix;

fname = fieldnames(SAT);

for i = 1 : length(SatID)
    [az,el,vis] = ...;
    Az_matrix(:,i) = az;
    El_matrix(:,i) = el;
    Vis_matrix(:,i) = vis;
end

El_matrix(El_matrix < 0) = missing;


%% SkyPlot 

figure();
sp = ...;
set(sp,'LabelFontSize',20);
for idx = 1:size(Az_matrix, 1)
    ...;
    drawnow limitrate
    pause(0.05);
end


%% Number of visible satellites

nSatView = ...;

plot(t/3600,nSatView);
xlabel('time [hr]'); ylabel('n. of visible s/c');

