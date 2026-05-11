function [sp3, Greg_time] = read_sp3_multiconstellation( filename , constID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign a file ID and open the given header file.
fid=fopen(filename);

% If the file does not exist, display warning message
if fid == -1
    display('Error!  File does not exist.  Computer will self-detonate');
else
    
    % Go through the header 
    F = 1;
    i = 1;
    while F(1) ~= '*'
        current_line = fgetl(fid);
        % Store the number of satellites in the SP3 file    
        if i == 3
            current_line = current_line(2:length(current_line));
            F = sscanf(current_line,'%u');
            no_sat = F(1);
        else
            F = sscanf(current_line,'%c');
        end
        i = i + 1;
        %pause()
    end
    
    % Begin going through times and observations
    end_of_file = 0;
    i = 0; j = 1;
    while end_of_file ~= 1        
        current_line = current_line(2:length(current_line));
        F = sscanf(current_line,'%f');
        % Load GPS Gregorian time into variables
        Y = F(1);
        M = F(2);
        D = F(3);
        H = F(4);
        min = F(5);
        sec = F(6);
        Greg_time(j,:) = [Y M D H min sec];
        
        % Convert GPS Gregorian time to GPS week and GPS TOW
        [GPS_wk, GPS_TOW] = GPSweek(Y,M,D,H,min,sec);
         
        % Store satellite PRN and appropriate observations
        for n = 1:no_sat
            
            % Go to the next line
            current_line = fgetl(fid);
            F = sscanf(current_line(3:length(current_line)),'%f');
            
            % Save PRN, positions, and clock error
            constellation = current_line(2);
            PRN = F(1); x = F(2); y = F(3); z = F(4); clk_err = F(5);
            
            % Create observation vector
            if constellation == 'G'
                sp3_obs_GPS(i+n,:) = table(datetime(Greg_time(j,:)),GPS_wk,GPS_TOW,PRN,x,y,z,clk_err);
            elseif constellation == 'E'
                sp3_obs_GALILEO(i+n,:) = table(datetime(Greg_time(j,:)),GPS_wk,GPS_TOW,PRN,x,y,z,clk_err);
            elseif constellation == 'R'
                sp3_obs_GLONASS(i+n,:) = table(datetime(Greg_time(j,:)),GPS_wk,GPS_TOW,PRN,x,y,z,clk_err);
            elseif constellation == 'C'
                sp3_obs_BEIDOU(i+n,:) = table(datetime(Greg_time(j,:)),GPS_wk,GPS_TOW,PRN,x,y,z,clk_err);
            elseif constellation == 'J'
                sp3_obs_QZSS(i+n,:) = table(datetime(Greg_time(j,:)),GPS_wk,GPS_TOW,PRN,x,y,z,clk_err);
            end
            n = n + 1;
        end
        
        % Go to next line - check to see if it is the end of file
        current_line = fgetl(fid);
        if strfind(current_line,'EOF')
            end_of_file = 1;
        end
        
        i = i + n - 1;
        j = j + 1;
    end         
end
% Data from SP3 file in ECEF reference frame
if exist('sp3_obs_GPS','var') && constID == "GPS"
    sp3.data = sp3_obs_GPS(any(table2array(sp3_obs_GPS(:,2:end)),2),:);
end
if exist('sp3_obs_GALILEO','var') && constID == "GALILEO"
    sp3.data = sp3_obs_GALILEO(any(table2array(sp3_obs_GALILEO(:,2:end)),2),:); 
end
if exist('sp3_obs_GLONASS','var') && constID == "GLONASS"
    sp3.data = sp3_obs_GLONASS(any(table2array(sp3_obs_GLONASS(:,2:end)),2),:);
end
if exist('sp3_obs_BEIDOU','var') && constID == "BEIDOU"
    sp3.data = sp3_obs_BEIDOU(any(table2array(sp3_obs_BEIDOU(:,2:end)),2),:);
end
if exist('sp3_obs_QZSS','var') && constID == "QZSS"
    sp3.data = sp3_obs_QZSS(any(table2array(sp3_obs_QZSS(:,2:end)),2),:);
end


% SP3 data columns
sp3.col.Date = 1;       % Date
sp3.col.WEEK = 2;       % week
sp3.col.TOW = 3;        % time of week in seconds
sp3.col.PRN = 4;        % Satellite ID
sp3.col.X = 5;          % X coordinate [km]
sp3.col.Y = 6;          % Y coordinate [km]
sp3.col.Z = 7;          % Z coordinate [km]
sp3.col.B = 8;          % clock error [s]



