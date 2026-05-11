"""Helpers for working with SP3 precise orbit data loaded via gnss_lib_py."""
import numpy as np


def sp3_get_sc_pos(sp3, sat_ids):
    """Extract ECEF positions for the requested PRNs.

    Parameters
    ----------
    sp3     : gnss_lib_py NavData/Sp3   (already filtered by gnss_id if multi-constellation)
    sat_ids : list[int]                 PRN numbers

    Returns
    -------
    sat       : dict[int, np.ndarray]   PRN -> (N, 3) ECEF positions [km],
                                        ordered by SP3 epoch
    available : list[int]               PRNs that actually had data
    """
    sat, available = {}, []
    for prn in sat_ids:
        record = sp3.where("sv_id", prn)
        if record is None or len(record) == 0:
            print(f"[warn] PRN {prn} not in SP3 — skipping")
            continue
        # gnss_lib_py stores SP3 positions in metres → convert to km
        sat[prn] = np.column_stack([record["x_sv_m"],
                                    record["y_sv_m"],
                                    record["z_sv_m"]]) / 1000.0
        available.append(prn)
    return sat, available