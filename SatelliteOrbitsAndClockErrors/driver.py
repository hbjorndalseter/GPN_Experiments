"""
driver.py — top-level pipeline.

Loads broadcast (RINEX nav) + precise (SP3) products, propagates positions
and clock errors for a chosen set of PRNs onto a common time grid, and
returns two dicts (broadcast, precise) ready for downstream analysis.
"""
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
import numpy as np
import gnss_lib_py as glp

from time_utils.gps_time          import greg2gps
from brdc_utils.get_sc_data       import get_sc_data
from sp3_utils.sp3_helpers        import sp3_get_sc_pos
from orbit_computation.gps        import gps_coordinates, MU as MU_GPS
from orbit_computation.glonass    import glonass_coordinates
from orbit_computation.sp3_interp import interp_precise_orbits


# ===================== CONFIG =====================
RINEX_FILE = "rinex_obs_files/BRDC00IGS_R_20240310000_01D_MN.rnx"
SP3_FILE   = "sp3_files/ESA0MGNFIN_20240310000_01D_05M_ORB.SP3"

CONST_ID      = "GLONASS"                    # "GPS" / "GALILEO" / "BEIDOU" / "GLONASS"
SAT_IDS       = [1, 15, 30]
DATE_UTC      = [2024, 1, 31, 0, 0, 0] 
DURATION_HRS  = 24.0
TIME_STEP_SEC = 60.0

# Simple conversion
GNSS_ID = {"GPS": "gps", "GALILEO": "galileo",
           "BEIDOU": "beidou", "GLONASS": "glonass"}
F_REL = -4.442807633e-10
C     = 2.99792458e8


# ===================== HELPERS =====================
def _pick_ephem(records, t_mid):
    """Index of the ephemeris record whose absolute time is closest to t_mid (SOW)."""
    # Robustly derive the record's monotonic Second of Week from absolute gps_millis
    rec_sow = (np.asarray(records["gps_millis"]) / 1000.0) - (gps_week * 604800)
    return int(np.argmin(np.abs(rec_sow - t_mid)))


def _solve_kepler(M, e, n_iter=15):
    E = M.copy() if np.ndim(M) else float(M)
    for _ in range(n_iter):
        E = M + e * np.sin(E)
    return E


# ===================== 1. TIME GRID =====================
gps_week, t_start_sow, _, _ = greg2gps(DATE_UTC)
n_pts        = int(DURATION_HRS * 3600.0 / TIME_STEP_SEC) + 1
t_target_sow = t_start_sow + np.arange(n_pts) * TIME_STEP_SEC
t_mid_sow    = t_start_sow + DURATION_HRS * 1800.0       # mid-window pick


# ===================== 2. LOAD FILES =====================
print(f"Loading {RINEX_FILE} ...")
ALL_SATS = ([f"G{i:02d}" for i in range(1, 33)] +    # GPS
            [f"E{i:02d}" for i in range(1, 37)] +    # Galileo
            [f"R{i:02d}" for i in range(1, 25)] +    # GLONASS
            [f"C{i:02d}" for i in range(1, 47)])     # BeiDou

rinex_nav = glp.RinexNav(RINEX_FILE, satellites=ALL_SATS)
print(f"Loading {SP3_FILE} ...")
sp3 = glp.Sp3(SP3_FILE)

const_lc = GNSS_ID[CONST_ID]
sp3_c    = sp3.where("gnss_id", const_lc)


# ===================== 3. BROADCAST PROPAGATION =====================
broadcast = {}
sat_records, _ = get_sc_data(rinex_nav, SAT_IDS, CONST_ID)

for prn, recs in sat_records.items():
    n_targets = len(t_target_sow)
    pos_arr   = np.zeros((n_targets, 3))
    clk_arr   = np.zeros(n_targets)

    for i, t in enumerate(t_target_sow):
        idx = _pick_ephem(recs, t)
        row = lambda k: float(np.asarray(recs[k])[idx])

        if CONST_ID == "GLONASS":
            # 1. Coordinate Time Alignment
            # GLONASS uses Moscow Time (UTC+3). Leap seconds = 18.
            t_moscow_target = t - 18.0 + 10800.0

            # 2. Lookup the correct Ephemeris for THIS specific time 't'
            # We must redefine 'idx' and 'row' inside the loop for every 't'
            rec_times = np.asarray(recs["gps_millis"]) / 1000.0
            rec_sow   = rec_times - (gps_week * 604800)
            idx = int(np.argmin(np.abs(rec_sow - t)))
            
            row = lambda k: float(np.asarray(recs[k])[idx])
            toe_raw = row("t_oe") 

            # 3. Calculate integration delta (dt) and handle day rollover
            dt = t_moscow_target - toe_raw
            if dt > 43200:  dt -= 86400
            elif dt < -43200: dt += 86400

            # 4. Prepare Ephemeris Dictionary
            ephem = {
                "PositionX": row("X"), "PositionY": row("Y"), "PositionZ": row("Z"),
                "VelocityX": row("dX"), "VelocityY": row("dY"), "VelocityZ": row("dZ"),
                "AccelerationX": row("dX2"), "AccelerationY": row("dY2"), "AccelerationZ": row("dZ2"),
                "t_oe": toe_raw, 
            }
            
            # 5. Call your Runge-Kutta function
            # We calculate position at (toe_raw + dt)
            try:
                x_km, y_km, z_km = glonass_coordinates(ephem, toe_raw + dt)
                
                # Convert km to m for comparison with SP3
                pos_arr[i] = [np.squeeze(x_km) * 1000.0, 
                              np.squeeze(y_km) * 1000.0, 
                              np.squeeze(z_km) * 1000.0]
                
            except Exception as e:
                pos_arr[i] = [np.nan, np.nan, np.nan]
            
            # 6. GLONASS Clock Correction
            tau_n = -row("SVclockBias")  # Note the negative sign common in RINEX
            gamma_n = row("SVrelFreqBias")
            clk_arr[i] = tau_n + gamma_n * dt

        else:  # GPS / Galileo / BeiDou
            e, sqrt_a  = row("e"), row("sqrtA")
            t_oe, t_oc = row("t_oe"), row("t_oc")
            delta_n    = row("deltaN")
            M0         = row("M_0")

            # Explicitly align evaluation time to BeiDou Time (BDT = GPST - 14s)
            t_eval = t - 14.0 if CONST_ID == "BEIDOU" else t

            # Force t_eval explicitly into the propagation function
            X, Y, Z = gps_coordinates(
                t_eval,  
                row("C_rs"), delta_n, M0,
                row("C_uc"), e, row("C_us"), sqrt_a, t_oe,
                row("C_ic"), row("Omega_0"), row("C_is"),
                row("i_0"), row("C_rc"), row("omega"),
                row("OmegaDot"), row("IDOT"),
            )
            pos_arr[i] = [np.squeeze(X), np.squeeze(Y), np.squeeze(Z)]

            a   = sqrt_a**2
            n   = np.sqrt(MU_GPS / a**3) + delta_n
            Mk  = M0 + n * (t_eval - t_oe)
            Ek  = _solve_kepler(Mk, e)
            dt_clk  = t_eval - t_oc
            clk_arr[i] = (row("SVclockBias")
                          + row("SVclockDrift")     * dt_clk
                          + row("SVclockDriftRate") * dt_clk * dt_clk
                          + F_REL * e * sqrt_a * np.sin(Ek))

    broadcast[prn] = {"pos_m": pos_arr, "clock_s": clk_arr, "time_sow": t_target_sow.copy()}


# ===================== 4. PRECISE ORBIT + CLOCK =====================
precise = {}
sat_sp3_pos, avail_sp3 = sp3_get_sc_pos(sp3_c, SAT_IDS)

for prn in avail_sp3:
    recs = sp3_c.where("sv_id", prn)
    
    sp3_abs_sec = np.asarray(recs["gps_millis"]) / 1000.0
    
    base_gps_week = int(sp3_abs_sec[len(sp3_abs_sec) // 2] // 604800)
    
    sp3_sow_monotonic = sp3_abs_sec - (base_gps_week * 604800)

    # Interpolate precise positions safely using aligned, monotonic SOW grids
    pos = interp_precise_orbits(sat_sp3_pos[prn], sp3_sow_monotonic, t_target_sow)
    
    # Robustly check if precise clock bias exists and interpolate
    if "b_sv_m" in recs.rows:
        clk_s = np.asarray(recs["b_sv_m"]) / C
        clk   = np.interp(t_target_sow, sp3_sow_monotonic, clk_s)
    else:
        clk   = np.zeros_like(t_target_sow)

    precise[prn] = {"pos_m": pos, "clock_s": clk, "time_sow": t_target_sow.copy()}


# ===================== 5. QUICK SUMMARY =====================
print(f"Satellite System: {CONST_ID}")
print(f"\nBroadcast PRNs: {sorted(broadcast)}    Precise PRNs: {sorted(precise)}")
common = sorted(set(broadcast) & set(precise))
for prn in common:
    diff = np.linalg.norm(broadcast[prn]["pos_m"] - precise[prn]["pos_m"], axis=1)
    print(f"  PRN {prn:2d}: |broadcast − precise|  "
          f"mean = {diff.mean():7.2f} m   max = {diff.max():7.2f} m")