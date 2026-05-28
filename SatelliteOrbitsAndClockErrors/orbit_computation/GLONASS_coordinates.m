function [x, y, z] = GLONASS_coordinates(ephem, interp_time)
% Convert broadcast values from km to meters for numerical integration
pos_0 = [ephem.PositionX, ephem.PositionY, ephem.PositionZ]' * 1e3;   % m
vel_0 = [ephem.VelocityX, ephem.VelocityY, ephem.VelocityZ]' * 1e3;   % m/s
y0    = [pos_0(1); vel_0(1); pos_0(2); vel_0(2); pos_0(3); vel_0(3)];
acc_p = [ephem.AccelerationX, ephem.AccelerationY, ephem.AccelerationZ]' * 1e3;  % m/s²

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-9);
t_query = interp_time(:);
t_pos   = sort(t_query(t_query >= 0), 'ascend');
t_neg   = sort(t_query(t_query <  0), 'ascend');

x_all = zeros(length(t_query), 1);
y_all = zeros(length(t_query), 1);
z_all = zeros(length(t_query), 1);

% Forward integration
if ~isempty(t_pos)
    [t_f, y_f] = ode45(@(t,yv) GLONASS_eq(yv, acc_p), [0, max(t_pos)], y0, opts);
    fprintf('ode45 used %d internal steps for span [0, %.1f] s\n', length(t_f), max(t_pos));
    [t_f, y_f] = ode45(@(t,yv) GLONASS_eq(yv, acc_p), [0; t_pos], y0, opts);
    x_all(t_query >= 0) = interp1(t_f, y_f(:,1), t_pos);
    y_all(t_query >= 0) = interp1(t_f, y_f(:,3), t_pos);
    z_all(t_query >= 0) = interp1(t_f, y_f(:,5), t_pos);
end

% Backward integration
if ~isempty(t_neg)
    [t_b, y_b] = ode45(@(t,yv) GLONASS_eq(yv, acc_p), [0; min(t_neg)], y0, opts);
    [t_b_s, si] = sort(t_b, 'ascend');
    y_b_s = y_b(si, :);
    x_all(t_query < 0) = interp1(t_b_s, y_b_s(:,1), t_neg);
    y_all(t_query < 0) = interp1(t_b_s, y_b_s(:,3), t_neg);
    z_all(t_query < 0) = interp1(t_b_s, y_b_s(:,5), t_neg);
end

% Convert back to km for output
x = x_all / 1e3;
y = y_all / 1e3;
z = z_all / 1e3;
end