%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute GPS satellite ECEF position from broadcast ephemeris
% following the IS-GPS-200 algorithm.
%
% Inputs (all SI units as broadcast):
%   t            - GPS time of signal transmission [s]
%   crs          - Amplitude of sin harmonic correction to orbit radius [m]
%   delta_n      - Mean motion difference from computed value [rad/s]
%   M0           - Mean anomaly at reference time [rad]
%   cuc          - Amplitude of cos harmonic correction to arg. of latitude [rad]
%   e            - Eccentricity [-]
%   cus          - Amplitude of sin harmonic correction to arg. of latitude [rad]
%   sqrt_a       - Square root of semi-major axis [sqrt(m)]
%   toe          - Reference time of ephemeris [s]
%   cic          - Amplitude of cos harmonic correction to inclination [rad]
%   nodo         - Longitude of ascending node at weekly epoch [rad]
%   cis          - Amplitude of sin harmonic correction to inclination [rad]
%   i0           - Inclination angle at reference time [rad]
%   crc          - Amplitude of cos harmonic correction to orbit radius [m]
%   long_perigeo - Argument of perigee [rad]
%   nodo_dot     - Rate of right ascension [rad/s]
%   i_dot        - Rate of inclination angle [rad/s]
%
% Outputs:
%   Xk, Yk, Zk  - ECEF satellite coordinates [km]
%
% Coded by Anese Giovanni / completed for Workbook 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xk, Yk, Zk] = GPS_coordinates(t, crs, delta_n, M0, cuc, e, cus, ...
                                          sqrt_a, toe, cic, nodo, cis, i0, ...
                                          crc, long_perigeo, nodo_dot, i_dot)

% ── Constants ──────────────────────────────────────────────────────────
mu      = 398600500000000;      % [m³/s²]  WGS-84 gravitational parameter
omega_e = 7.2921151467e-5;      % [rad/s]  WGS-84 Earth rotation rate

% ── Step 1: Semi-major axis ────────────────────────────────────────────
a = sqrt_a^2;                   % [m]

% ── Step 2: Computed mean motion ───────────────────────────────────────
n0 = sqrt(mu / a^3);            % [rad/s]

% ── Step 3: Time from ephemeris reference epoch ────────────────────────
tk = t - toe;
% Handle GPS week crossover
if tk >  302400, tk = tk - 604800; end
if tk < -302400, tk = tk + 604800; end

% ── Step 4: Corrected mean motion ─────────────────────────────────────
n = n0 + delta_n;               % [rad/s]

% ── Step 5: Mean anomaly ───────────────────────────────────────────────
Mk = M0 + n * tk;               % [rad]

% ── Step 6: Eccentric anomaly (Kepler's equation, iterative) ──────────
Ek = Mk;                        % initial guess
for iter = 1:50
    Ek_new = Mk + e * sin(Ek);
    if abs(Ek_new - Ek) < 1e-12
        break;
    end
    Ek = Ek_new;
end
Ek = Ek_new;

% ── Step 7: True anomaly ───────────────────────────────────────────────
sin_vk = sqrt(1 - e^2) * sin(Ek) / (1 - e * cos(Ek));
cos_vk = (cos(Ek) - e)           / (1 - e * cos(Ek));
vk = atan2(sin_vk, cos_vk);

% ── Step 8: Argument of latitude ───────────────────────────────────────
phi_k = vk + long_perigeo;

% ── Step 9: Second-order harmonic corrections ──────────────────────────
sin2phi = sin(2 * phi_k);
cos2phi = cos(2 * phi_k);

delta_uk = cus * sin2phi + cuc * cos2phi;   % arg. of latitude correction [rad]
delta_rk = crs * sin2phi + crc * cos2phi;   % radius correction           [m]
delta_ik = cis * sin2phi + cic * cos2phi;   % inclination correction      [rad]

% ── Step 10: Corrected orbital elements ────────────────────────────────
uk = phi_k + delta_uk;                      % corrected arg. of latitude  [rad]
rk = a * (1 - e * cos(Ek)) + delta_rk;     % corrected radius            [m]
ik = i0 + delta_ik + i_dot * tk;           % corrected inclination       [rad]

% ── Step 11: Position in orbital plane ─────────────────────────────────
xk_prime = rk * cos(uk);
yk_prime = rk * sin(uk);

% ── Step 12: Corrected longitude of ascending node ─────────────────────
Omega_k = nodo + (nodo_dot - omega_e) * tk - omega_e * toe;

% ── Step 13: ECEF coordinates [m] ──────────────────────────────────────
Xk = xk_prime * cos(Omega_k) - yk_prime * cos(ik) * sin(Omega_k);
Yk = xk_prime * sin(Omega_k) + yk_prime * cos(ik) * cos(Omega_k);
Zk = yk_prime * sin(ik);

% ── Convert to km (consistent with SP3 / rest of pipeline) ─────────────
Xk = Xk / 1e3;
Yk = Yk / 1e3;
Zk = Zk / 1e3;

end