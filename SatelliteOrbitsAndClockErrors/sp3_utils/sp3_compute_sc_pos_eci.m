%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code updates the sp3 structure with position of desired
% satellites in ECI reference frame
%
% Inputs:   sp3 - structure with info about sp3 data [GPS_week GPS_TOW PRN x y z]
%           Greg_time - time span vector [year,month,day,hour,minutem,second]
%           constID - ID of the constellation (es. GPS,GLONASS,ecc)
%
% Outputs:  sp3 - updated strucure with ECI positions
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sp3 = sp3_compute_sc_pos_eci(sp3,Greg_time, constID)

no_sat = height(sp3.data(:,1))/length(Greg_time); % number of satellites
sp3.(strcat(constID,'_ECI')) = sp3.data;   % preallocation
sp3.(strcat(constID,'_ECI')){:,5:7} = NaN(height(sp3.data(:,1)),3);

for i = 1:length(Greg_time)
    for j = 1 : no_sat
        % Extract ECEF position
        r_ecef = table2array(sp3.data(j+no_sat*(i-1),5:7));
        % Conversion from ECEF to ECI
        r_eci = ecef2eci(Greg_time(i,:),r_ecef);
        % Update structure with ECI position
        sp3.(strcat(constID,'_ECI')){j+no_sat*(i-1),5:7} = r_eci';
    end
end

end