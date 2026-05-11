% Exercise 1: GNSS 2D
clear, clc, close all
problem_type = 1; % 1 = without the receiver clock bias
example_number = 1; % 1 = good position of the ground stations; 2 = stations alligned
[x0_true, y0_true, GS, GS_n, p_meas, Dx, R_pos, ctau_true, R_ctau] = data_generation(problem_type,example_number);
%% main while loop
% loop parameters
tol = 0.1; % tolerance or thershold fro the while loop end
stop_check = 1; % evaluated distance from previous position
k_max = 15; % max number oof iteration
k = 0; % interation counter
% variables initialization
d_calc = zeros(GS_n,1);
p_calc = zeros(GS_n,1);
Positions = zeros(k_max,2);
Positions(1,:) = R_pos;
r_norm_pre = zeros(k_max,1); % norm of the pre fit residuals
r_norm_post = zeros(k_max,1);% norm of the post fit residuals
if problem_type == 1
    H = zeros(GS_n,2);
    Dx_tot = zeros(k_max,2);
end
while ((stop_check>tol)&&(k<k_max))
    k = k+1;
    % distances and pseudoranges
    for i = 1:GS_n
        d_calc(i) = sqrt((GS(i,1)-R_pos(1))^2+(GS(i,2)-R_pos(2))^2);
        p_calc(i) = d_calc(i);
    end
    % pre fit residuals
    Y_res = p_meas - p_calc;
    r_norm_pre(k) = Y_res'*Y_res; % squared value of the norm
    % geometry matrix H
    for i=1:GS_n
        H(i,1) = -(GS(i,1)-R_pos(1))/d_calc(i);
        H(i,2) = -(GS(i,2)-R_pos(2))/d_calc(i);
    end
    % evaluation of the solution
    AA = H'*H;
    bb = H'*Y_res;
    Dx_new = AA\bb;
    % stop condition evaluation
    stop_check = norm(Dx_new);
    % new receiver position
    Dx = Dx_new;
    R_pos = R_pos + Dx';
    Positions(k+1,:) = R_pos;
    Dx_tot(k,:) = Dx';
end

