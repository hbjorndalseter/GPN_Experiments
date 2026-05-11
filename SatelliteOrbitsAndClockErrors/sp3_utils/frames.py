"""ECEF <-> ECI rotation via Greenwich Mean Sidereal Time.

This is a single-rotation model (Earth spin only — no precession, nutation,
polar motion, or UT1-UTC). Good to ~arcsecond level for the day-scale
workbook orbits. For sub-cm-grade work, swap to astropy.coordinates.
"""
import numpy as np
from datetime import datetime, timezone


def _gmst_radians(dt_utc):
    """GMST at a UTC datetime, in radians (Vallado Alg. 15, ~1 arcsec accuracy)."""
    j2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
    if dt_utc.tzinfo is None:
        dt_utc = dt_utc.replace(tzinfo=timezone.utc)
    T = (dt_utc - j2000).total_seconds() / 86400.0 / 36525.0    # Julian centuries
    gmst_sec = (67310.54841
                + (876600.0 * 3600.0 + 8640184.812866) * T
                + 0.093104 * T**2
                - 6.2e-6 * T**3)
    return np.deg2rad((gmst_sec / 240.0) % 360.0)                # 240 sec per deg


def ecef2eci(greg_time, r_ecef):
    """Rotate ECEF position(s) to ECI.

    greg_time : 6-tuple/array [Y, M, D, H, min, sec]  or  datetime  (UTC)
    r_ecef    : (3,) or (N, 3) ECEF vector(s)
    Returns the same shape as r_ecef.
    """
    if isinstance(greg_time, datetime):
        dt = greg_time
    else:
        y, mo, d, h, mi, s = greg_time
        sec_int, sec_frac = divmod(float(s), 1.0)
        dt = datetime(int(y), int(mo), int(d), int(h), int(mi),
                      int(sec_int), int(sec_frac * 1e6), tzinfo=timezone.utc)

    theta = _gmst_radians(dt)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s, 0.0],
                  [s,  c, 0.0],
                  [0.0, 0.0, 1.0]])

    r = np.atleast_2d(r_ecef)
    out = r @ R.T                       # (N,3) @ (3,3) -> (N,3)
    return out[0] if np.ndim(r_ecef) == 1 else out


def sp3_compute_sc_pos_eci(ecef_positions, greg_times):
    """Convert an aligned series of SP3 ECEF positions to ECI.

    ecef_positions : (N, 3) [km]
    greg_times     : (N, 6) [Y, M, D, H, min, sec] UTC, one row per ECEF row
    """
    ecef_positions = np.asarray(ecef_positions, dtype=float)
    greg_times = np.atleast_2d(np.asarray(greg_times, dtype=float))
    out = np.empty_like(ecef_positions)
    for i in range(len(greg_times)):
        out[i] = ecef2eci(greg_times[i], ecef_positions[i])
    return out