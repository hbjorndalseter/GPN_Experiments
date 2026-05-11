"""GPS / Galileo / BeiDou Keplerian broadcast orbit propagation (ICD-GPS-200)."""
import numpy as np

MU      = 3.986005e14
OMEGA_E = 7.2921151467e-5


def gps_coordinates(t, crs, delta_n, M0, cuc, e, cus, sqrt_a, toe,
                    cic, Omega0, cis, i0, crc, omega, OmegaDot, IDOT):
    """Compute satellite ECEF position from broadcast Keplerian elements.

    Argument order and names mirror MATLAB GPS_coordinates(), with renames:
        nodo         -> Omega0       (RAAN at week start)
        long_perigeo -> omega        (argument of perigee)
        nodo_dot     -> OmegaDot     (rate of RAAN)
        i_dot        -> IDOT         (rate of inclination)

    t may be scalar or array. Returns Xk, Yk, Zk in metres.
    """
    t = np.asarray(t, dtype=float)
    a = sqrt_a**2
    n = np.sqrt(MU / a**3) + delta_n

    tk = t - toe
    tk = np.where(tk >  302400.0, tk - 604800.0, tk)
    tk = np.where(tk < -302400.0, tk + 604800.0, tk)

    Mk = M0 + n*tk
    # Solve Kepler's equation by fixed-point iteration; converges in <10 for GNSS e
    Ek = Mk.copy() if np.ndim(Mk) else float(Mk)
    for _ in range(15):
        Ek = Mk + e*np.sin(Ek)

    sinE, cosE = np.sin(Ek), np.cos(Ek)
    vk  = np.arctan2(np.sqrt(1.0 - e*e)*sinE, cosE - e)
    Phi = vk + omega

    s2, c2 = np.sin(2*Phi), np.cos(2*Phi)
    uk = Phi      + cus*s2 + cuc*c2
    rk = a*(1.0 - e*cosE) + crs*s2 + crc*c2
    ik = i0 + IDOT*tk     + cis*s2 + cic*c2

    xp, yp = rk*np.cos(uk), rk*np.sin(uk)
    Omk = Omega0 + (OmegaDot - OMEGA_E)*tk - OMEGA_E*toe

    Xk = xp*np.cos(Omk) - yp*np.cos(ik)*np.sin(Omk)
    Yk = xp*np.sin(Omk) + yp*np.cos(ik)*np.cos(Omk)
    Zk = yp*np.sin(ik)
    return Xk, Yk, Zk