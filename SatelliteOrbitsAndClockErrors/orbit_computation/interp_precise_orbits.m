function [orbit] = interp_precise_orbits(pos, time, interp_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolates satellite positions from SP3 file data using
% Lagrange polynomial interpolation (standard IGS method, 10 points).
%
% Inputs:
%   pos         - Nx3 matrix of satellite ECEF coordinates [km]
%   time        - Nx1 time vector corresponding to SP3 epochs [s]
%   interp_time - Mx1 vector of query times [hours]
%
% Outputs:
%   orbit       - Mx3 matrix of interpolated ECEF coordinates [km]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interp_time_s = interp_time(:) * 3600;   % hours → seconds
n_query       = length(interp_time_s);
orbit         = zeros(n_query, 3);
npts          = 10;                       % Lagrange window size (IGS standard)

for i = 1:n_query
    t_q = interp_time_s(i);

    [~, idx] = min(abs(time - t_q));

    % Build symmetric window, clamped to valid range
    half    = floor(npts / 2);
    i_start = max(1, idx - half);
    i_end   = i_start + npts - 1;
    if i_end > length(time)
        i_end   = length(time);
        i_start = max(1, i_end - npts + 1);
    end

    t_win = time(i_start:i_end);

    for c = 1:3
        orbit(i, c) = lagrange_eval(t_win, pos(i_start:i_end, c), t_q);
    end
end

end

% ── Lagrange polynomial evaluation ────────────────────────────────────
function y = lagrange_eval(t, p, t_q)
n = length(t);
y = 0;
for i = 1:n
    Li = 1;
    for j = 1:n
        if j ~= i
            Li = Li * (t_q - t(j)) / (t(i) - t(j));
        end
    end
    y = y + p(i) * Li;
end
end