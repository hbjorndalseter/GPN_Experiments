%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by Andrea Valmorbida and Anese Giovanni
% DII and CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% workspace setting

close all; clear; clc;

addpath(genpath('.\'));


%% Initial settings

% file available at https://cddis.nasa.gov/archive/gnss/data/daily/
filename = 'BRDC00IGS_R_20240310000_01D_MN.rnx';

% Selected constellation
constID = "GPS";

% Satellite ID (choose satellite ID to study)
SatID = [1,15,30];

% Time to interpolate
vector_time = linspace(0,3,1000);  % [hours]
toe_start = 259200;


%%  Read data from RINEX OBSERVATION FILE

% Read broadcast file
brdc = ...;
info = ...;

% Ephemeris of the selected constellation
brdcEphem = brdc.(constID);


%% Study selected satellites

[SAT,SatID] = brdc_get_sc_data(...);


%% Orbit computation 

orbit = zeros(length(vector_time),3);
fname = fieldnames(SAT);

if constID == "GPS"
    for i = 1 : length(SatID)
        t = toe_start +interp_time*3600;
        % Compute the positions in the interpolation time points
        for j = 1 : length(interp_time)
            % Ephemeris data
            ...
            % get coordinates
            [x,y,z] = GPS_coordinates(...);
            orbit(j,:) = (...) ;
        end
        % Update SAT structure
        SAT.(string([char(fname(i)),'_brdc']))=orbit;
    end

% elseif constID == "GLONASS"
%     % Get info from RINEX file
%     info = ...;
%     % Leap seconds
%     Leap_Seconds = ...; % [s]
%     for i = 1 : length(SatID)
%         % Time
%         t = interp_time*3600-Leap_Seconds;
%         % Epehemeris data
%         ...
%         % get coordinates
%         [x,y,z] = GLONASS_coordinates(ephem(1,:), t);
%         orbit = ...;    % [km]   
%         % Update SAT structure
%         ...
%     end
end


%% Visualize orbits

...

