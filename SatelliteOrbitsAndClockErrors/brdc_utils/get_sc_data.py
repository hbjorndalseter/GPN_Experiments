import gnss_lib_py as glp

CONST_MAP = {"GPS": "gps", "GALILEO": "galileo", "GLONASS": "glonass",
             "BEIDOU": "beidou", "QZSS": "qzss"}

def get_sc_data(rinex_nav, sat_ids, const_id):
    """
    Extract per-satellite broadcast ephemeris records.

    Parameters
    ----------
    rinex_nav : glp.RinexNav        parsed RINEX navigation file
    sat_ids   : list[int]           PRN numbers you want
    const_id  : str                 "GPS", "GALILEO", "GLONASS", "BEIDOU", "QZSS"

    Returns
    -------
    sat       : dict[int, NavData]  PRN -> ephemeris record(s) for that satellite
    available : list[int]           PRNs that actually had data
    """
    const_subset = rinex_nav.where("gnss_id", CONST_MAP[const_id.upper()])

    sat = {}
    available = []
    for prn in sat_ids:
        record = const_subset.where("sv_id", prn)
        if record is not None and len(record) > 0:
            sat[prn] = record
            available.append(prn)
        else:
            print(f"[warn] {const_id} PRN {prn} not in broadcast file — skipping")
    return sat, available