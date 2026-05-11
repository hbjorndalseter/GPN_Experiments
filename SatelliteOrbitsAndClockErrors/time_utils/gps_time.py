"""GPS <-> Gregorian time conversions.

GPS epoch: 1980-01-06 00:00:00 UTC.

Caveat (matches the MATLAB originals): leap seconds are ignored — GPS time
and UTC are treated as identical. Fine for workbook exercises; if you ever
need leap-second-correct conversion use gnss_lib_py's helpers or astropy.time.
"""
from datetime import datetime, timedelta
import numpy as np

GPS_EPOCH    = datetime(1980, 1, 6, 0, 0, 0)
SEC_PER_WEEK = 7 * 86400


def greg2gps(time):
    """Calendar date(s) -> (gps_week, tow, doy, dow).

    Parameters
    ----------
    time : (6,) or (N, 6) array-like
        Columns: [year, month, day, hour, minute, second].

    Returns
    -------
    Scalars (int, float, int, int) for single input,
    NumPy arrays for batch input.
        gps_week : GPS week number (full, not mod 1024)
        tow      : seconds of week
        doy      : day of year (1–366)
        dow      : day of week (0=Sunday … 6=Saturday)
    """
    arr = np.atleast_2d(np.asarray(time, dtype=float))
    weeks = np.empty(len(arr), dtype=int)
    tows  = np.empty(len(arr), dtype=float)
    doys  = np.empty(len(arr), dtype=int)
    dows  = np.empty(len(arr), dtype=int)

    for k, (y, mo, d, h, mi, s) in enumerate(arr):
        dt   = datetime(int(y), int(mo), int(d)) + timedelta(
                   hours=h, minutes=mi, seconds=s)
        secs = (dt - GPS_EPOCH).total_seconds()
        weeks[k] = int(secs // SEC_PER_WEEK)
        tows[k]  = secs - weeks[k] * SEC_PER_WEEK
        doys[k]  = dt.timetuple().tm_yday
        dows[k]  = int(tows[k] // 86400)

    if arr.shape[0] == 1:
        return int(weeks[0]), float(tows[0]), int(doys[0]), int(dows[0])
    return weeks, tows, doys, dows


def GPSweek(Y, M, D, H=0, minute=0, sec=0):
    """Single date -> (gps_week, tow). Wrapper matching the MATLAB signature."""
    w, t, _, _ = greg2gps([Y, M, D, H, minute, sec])
    return w, t


def gps2greg(gps_week, gps_second):
    """(gps_week, gps_second) -> datetime. Accepts scalars or arrays."""
    gw = np.atleast_1d(gps_week).astype(int)
    gs = np.atleast_1d(gps_second).astype(float)
    out = np.array([GPS_EPOCH + timedelta(weeks=int(w), seconds=float(s))
                    for w, s in zip(gw, gs)])
    return out[0] if out.size == 1 else out