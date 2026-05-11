import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import gnss_lib_py as glp

# Load precise orbit files utilizing the Extended Standard Product- 3 (SP3c) format
filename = 'ESA0MGNFIN_20240310000_01D_05M_ORB.SP3'
sp3_data = glp.Sp3(filename)

# Filter ECEF positions ('x_sv_m', 'y_sv_m', 'z_sv_m') and timing ('gps_millis') for a target PRN
prn_1_sp3 = sp3_data.where("gnss_id", "gps")
if sp3_data is not None:
    prn_1_sp3 = sp3_data.where("sv_id", 1)

# Extract the true precise anchor points for cubic spline interpolation
# (Converting native milliseconds to hours, and native meters to kilometers)
if prn_1_sp3 is not None:
    sp3_time_grid = prn_1_sp3["gps_millis"] / 3600000.0 
    sp3_x = prn_1_sp3["x_sv_m"] / 1000.0  
    sp3_y = prn_1_sp3["y_sv_m"] / 1000.0  
    sp3_z = prn_1_sp3["z_sv_m"] / 1000.0  

# --- Output Directory Setup ---
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# --- Configuration ---
filename = 'ESA0MGNFIN_20240310000_01D_05M_ORB.SP3'
SatID = [1, 15, 30]
interp_time = np.linspace(0, 3, 1000)  # Dense interpolation target [hours]

# Mocking sparse SP3 extraction (1 epoch every 15 minutes = 0.25 hours)
sp3_time_grid = np.linspace(0, 3, 13)

SAT_PRECISE = {}

fig, ax = plt.subplots(figsize=(10, 6))
ax.set_title("SP3 Precise Orbit Interpolation Accuracy", fontsize=14, fontweight='bold')

for prn in SatID:
    # Generate mock precise anchor points exhibiting typical MEO radius (~26560 km)
    phase = prn * 0.5
    sp3_x = 26560.0 * np.cos(sp3_time_grid + phase)
    sp3_y = 26560.0 * np.sin(sp3_time_grid + phase)
    sp3_z = 15000.0 * np.sin(sp3_time_grid * 0.5)
    
    # Fit Cubic Splines independently for each spatial dimension
    cs_x = CubicSpline(sp3_time_grid, sp3_x)
    cs_y = CubicSpline(sp3_time_grid, sp3_y)
    cs_z = CubicSpline(sp3_time_grid, sp3_z)
    
    # Evaluate at highly dense interpolation targets
    dense_x = cs_x(interp_time)
    dense_y = cs_y(interp_time)
    dense_z = cs_z(interp_time)
    
    SAT_PRECISE[f"PRN_{prn}_sp3"] = np.column_stack((dense_x, dense_y, dense_z))
    
    # Plot dimension X interpolation validation
    p = ax.plot(interp_time, dense_x, linestyle='-', linewidth=1.5, label=f"PRN {prn} (Interp)")
    ax.scatter(sp3_time_grid, sp3_x, color=p[0].get_color(), marker='o', s=30, zorder=5)

ax.set_xlabel("Time [hours]", fontsize=12)
ax.set_ylabel("ECEF X Coordinate [km]", fontsize=12)
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend()
plt.tight_layout()

save_path = os.path.join(output_dir, "sp3_precise_interpolation.png")
plt.savefig(save_path, dpi=300)
print(f"Precise interpolation plot successfully saved to: {save_path}")
plt.show()