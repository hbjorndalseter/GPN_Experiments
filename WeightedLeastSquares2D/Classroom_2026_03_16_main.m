% Exercise 1: GNSS 2D
clear %clc, close all
problem_type = 2; % 1 = without the receiver clock bias; 2 = with receiver clock bias
example_number = 1; % 1 = good position of the ground stations; 2 = stations alligned
% different type of weights
weighted = 1; % 0 = no weights; 1 = weighted, uncorrelated measuremens and identical variance;
% 2 = weighted, uncorrelated measurements
[x0_true, y0_true, GS, GS_n, p_meas, Dx, R_pos, ctau_true, R_ctau] = data_generation(problem_type,example_number);
%% main while loop
% loop parameters
tol = 0.1; % tolerance or thershold fro the while loop end
stop_check = 1; % evaluated distance from previous position
k_max = 15; % max number of iterations
k = 0; % interation counter
% variables initialization
d_calc = zeros(GS_n,1);
p_calc = zeros(GS_n,1);
Positions = zeros(k_max,2);
Positions(1,:) = R_pos;
r_norm_pre = zeros(k_max,1); % norm of the pre fit residuals
r_norm_post = zeros(k_max,1);% norm of the post fit residuals
mean_Y = zeros(k_max,1);
var_Y = zeros(k_max,1);
if problem_type == 1
    H = zeros(GS_n,2);
    Dx_tot = zeros(k_max,2);
    Cxx_tot = zeros(2,2,k_max);
    Q_tot = zeros(2,2,k_max);
elseif problem_type == 2
    H = zeros(GS_n,3);
    estimated_ctau = zeros(k_max,1);
    estimated_ctau(1) = R_ctau;
    Dx_tot = zeros(k_max,3);
    Cxx_tot = zeros(3,3,k_max);
    Q_tot = zeros(3,3,k_max);
end
while ((stop_check>tol)&&(k<k_max))
    k = k+1;
    % distances and pseudoranges
    for i = 1:GS_n
        d_calc(i) = sqrt((GS(i,1)-R_pos(1))^2+(GS(i,2)-R_pos(2))^2);
        p_calc(i) = d_calc(i) + R_ctau;
    end
    % pre fit residuals
    Y_res = p_meas - p_calc;
    r_norm_pre(k) = Y_res'*Y_res; % squared value of the norm
    mean_Y(k) = mean(Y_res);
    var_Y(k) = 1/(GS_n-1)*(Y_res'*Y_res);
    % geometry matrix H
    for i=1:GS_n
        H(i,1) = -(GS(i,1)-R_pos(1))/d_calc(i);
        H(i,2) = -(GS(i,2)-R_pos(2))/d_calc(i);
        if problem_type == 2
            H(i,3) = 1;
        end
    end
    % Weights and Cpp
    if weighted == 0
        W = eye(GS_n);
        % the variance of Y is approximated by:
        Cpp = var_Y(k)*eye(GS_n);
    elseif weighted == 1
        Cpp = var_Y(k)*eye(GS_n);
        W = Cpp\eye(GS_n);%inv(Cpp)
    elseif weighted == 2
        Cpp = diag(Y_res.^2); % only a bad approximation using only one measurement for ech component
        W = Cpp\eye(GS_n);%inv(Cpp)
    end
    % evaluation of the solution
    AA = H'*W*H;
    bb = H'*W*Y_res;
    Dx_new = AA\bb;
    % stop condition evaluation
    stop_check = norm(Dx_new(1:2));
    % new receiver position
    Dx = Dx_new;
    R_pos = R_pos + Dx(1:2)';
    Positions(k+1,:) = R_pos;
    Dx_tot(k,:) = Dx';
    % clock update
    if problem_type == 2
        R_ctau = R_ctau + Dx(3);
        estimated_ctau(k+1) = R_ctau;
    end
    % Evaluation of post-fit residuals
    r_norm_post(k) = (Y_res-H*Dx)'*(Y_res-H*Dx)% squared value of the norm
    % matrix Q
    Q = (H'*H)\eye(size(H,2));
    Q_tot(:,:,k) = Q;
    % Covariance matrix of the solution
    if weighted == 0
        Cxx = Q*H'*Cpp*H*Q;
    else
        Cxx = AA\eye(size(H,2));% inv(H'*W*H)
    end
    Cxx_tot(:,:,k) = Cxx;
end
% Simplify vectors and matrices
Positions = Positions(1:k+1,:);
r_norm_pre = r_norm_pre(1:k);
r_norm_post = r_norm_post(1:k);
mean_Y = mean_Y(1:k);
var_Y = var_Y(1:k);
Dx_tot = Dx_tot(1:k,:);
Q_tot = Q_tot(:,:,1:k);  
Cxx_tot = Cxx_tot(:,:,1:k);
if problem_type == 2
    estimated_ctau = estimated_ctau(1:k+1);
end
% Charts
figure
plot(Positions(:,1),Positions(:,2),'k+--')
xlabel('X [length units]')
ylabel('Y [length units]')
title('Positions')
axis equal
hold on
grid on
plot(Positions(1,1),Positions(1,2),'mo--','LineWidth',2)
plot(Positions(k+1,1),Positions(k+1,2),'bx--','LineWidth',2)
plot(GS(:,1),GS(:,2),'r*')
plot(x0_true,y0_true,'g+','LineWidth',2)
legend('Positions','Initial','Final','Stations','true')
if problem_type == 2
    figure
    plot(0:k,estimated_ctau,'r--x')
    hold on
    grid on
    plot(0:k,ctau_true*ones(size(estimated_ctau)),'b')
    xlabel('Iterations')
    ylabel('length units')
    title('Estimated ctau')
    legend('Estimated','treu')
end
figure
plot(Positions(:,1),Positions(:,2),'k+--')
title('Positions with uncertainties')
axis equal
hold on
grid on
plot(Positions(1,1),Positions(1,2),'mo--','LineWidth',2)
plot(Positions(k+1,1),Positions(k+1,2),'bx--','LineWidth',2)
plot(x0_true,y0_true,'g+','LineWidth',2)
LevConf = 0.9889; % K = 3 for a 2D problem
for i = 1:k
    error_ellipse(Cxx_tot(1:2,1:2,i),Positions(i+1,:)','conf',LevConf)
end
% Chart of residuals
figure
plot(0:k-1,r_norm_pre,'r--o')
hold on
grid on
plot(1:k,r_norm_post,'b--x')
title('Prefit and postfit residuals')
xlabel('Iterations')
ylabel('[cm^2]')
legend('Prefit','Postfit')


