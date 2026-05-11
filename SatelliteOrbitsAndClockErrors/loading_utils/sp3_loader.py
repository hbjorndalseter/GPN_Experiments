# loading_utils/sp3_loader.py
import gnss_lib_py as glp

def load_sp3(filename, const_id="GPS"):
    """Mimics MATLAB read_sp3_multiconstellation interface using gnss_lib_py."""
    const_map = {"GPS": "gps", "GALILEO": "galileo", "GLONASS": "glonass",
                 "BEIDOU": "beidou", "QZSS": "qzss"}
    sp3 = glp.Sp3(filename).where("gnss_id", const_map[const_id.upper()])
    # Return a simple namespace matching the MATLAB columns:
    # gps_millis, sv_id, x_sv_m, y_sv_m, z_sv_m  (units: ms and m, not s and km)
    return sp3