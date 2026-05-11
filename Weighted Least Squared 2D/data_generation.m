function [x0_true, y0_true, GS, GS_n, p_meas, Dx, R_pos, ctau_true, R_ctau] = data_generation(problem_type,example_number)
% receiver true postion
x0_true = 4.4;
y0_true = 3.0;
% Initial position error
ex = -2;
ey = 1.5;
% initial receiver position
R_x0 = x0_true+ex;
R_y0 = y0_true+ey;
R_pos = [R_x0 R_y0];
switch example_number
    case 1  
        % positions of ground stations
        GS_x = [9.0 0.2 3.7 -9.2 -5.5 8.4 -6.0 10.4 10.5]';
        GS_y = [4.0 7.7 1.5 5.8 -2.8 -5.7 1.0 -0.8 8.0]';
        % measured distances (pseudoranges)
        p_meas = [4.7 6.3 1.7 13.9 11.5 9.6 10.6 7.2 7.9]';
    case 2
        % positions of ground stations
        GS_x = [9.8 6.2 6.7 0.6 3.3 -2.0 -2.8]';
        GS_y = [7.8 6.7 3.2 2.5 -2.7 -1.0 -6.0]';
        % measured distances (pseudoranges)
        p_meas = [7.2 4.1 2.3 3.8 5.8 7.6 11.5]';
end
if problem_type == 1
    Dx = [0 0]';
    ctau_true = 0;
    R_ctau = 0;
end
% matrix of ground station positions
GS = [GS_x GS_y];
% number of ground stations
GS_n = length(GS_x);