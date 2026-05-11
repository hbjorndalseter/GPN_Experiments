"""GLONASS broadcast orbit propagation (PZ-90 ECEF frame).

Units: km, km/s, km/s^2.
Matches Workbook 2 task spec and the native units of GLONASS RINEX records.
"""
import numpy as np
from scipy.integrate import solve_ivp

# Task-specified constants (Workbook 2)
MU      = 3.9860044e5     # km^3 / s^2
OMEGA_E = 0.7292115e-4    # rad / s
J2      = -1.08263e-3     # NB: negative — paired with +(3/2)*J2*... in glonass_eq
RE      = 6.378136e3      # km


def glonass_eq(t, y, acc_ls):
    """RHS of GLONASS equations of motion in PZ-90 ECEF.

    y      = [x, vx, y, vy, z, vz]   (km, km/s)
    acc_ls = (ax, ay, az)            lunisolar accelerations (km/s^2),
                                     constant over the broadcast validity window.
    """
    x, vx, yp, vy, z, vz = y
    r2 = x*x + yp*yp + z*z
    r  = np.sqrt(r2)
    r3, r5 = r2*r, r2*r2*r

    # j2_pref is negative because J2 itself is negative in this convention.
    # The leading sign in each acceleration line below stays "+" as in the workbook.
    j2_pref = 1.5 * J2 * MU * RE**2 / r5
    z_xy = 1.0 - 5.0*z*z/r2
    z_zz = 3.0 - 5.0*z*z/r2
    ax_ls, ay_ls, az_ls = acc_ls

    return [
        vx,
        -MU*x/r3  + j2_pref*x*z_xy  + OMEGA_E**2*x  + 2*OMEGA_E*vy + ax_ls,
        vy,
        -MU*yp/r3 + j2_pref*yp*z_xy + OMEGA_E**2*yp - 2*OMEGA_E*vx + ay_ls,
        vz,
        -MU*z/r3  + j2_pref*z*z_zz  + az_ls,
    ]


def glonass_coordinates(ephem, interp_time_s):
    """Propagate a GLONASS satellite from broadcast ephemeris.

    Parameters
    ----------
    ephem : dict-like (units: km, km/s, km/s^2)
        PositionX/Y/Z, VelocityX/Y/Z, AccelerationX/Y/Z, t_oe (seconds)
    interp_time_s : target times [s], same epoch convention as t_oe

    Returns
    -------
    x, y, z : ECEF coordinate arrays [km]
    """
    y0 = np.array([ephem["PositionX"], ephem["VelocityX"],
                   ephem["PositionY"], ephem["VelocityY"],
                   ephem["PositionZ"], ephem["VelocityZ"]], dtype=float)
    acc   = (ephem["AccelerationX"], ephem["AccelerationY"], ephem["AccelerationZ"])
    t_ref = ephem["t_oe"]
    t_eval = np.atleast_1d(np.asarray(interp_time_s, dtype=float))

    fwd, bwd = t_eval >= t_ref, t_eval < t_ref
    out = np.empty((6, t_eval.size))

    if fwd.any():
        sol = solve_ivp(lambda t, y: glonass_eq(t, y, acc),
                        (t_ref, t_eval[fwd].max()), y0,
                        t_eval=t_eval[fwd],
                        method="RK45", rtol=1e-12, atol=1e-9, max_step=60.0)
        out[:, fwd] = sol.y
    if bwd.any():
        sol = solve_ivp(lambda t, y: glonass_eq(t, y, acc),
                        (t_ref, t_eval[bwd].min()), y0,
                        t_eval=t_eval[bwd][::-1],
                        method="RK45", rtol=1e-12, atol=1e-9, max_step=60.0)
        out[:, bwd] = sol.y[:, ::-1]

    return out[0], out[2], out[4]