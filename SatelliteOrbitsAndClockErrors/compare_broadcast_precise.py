import os
import numpy as np
import matplotlib.pyplot as plt

# --- Output Directory Setup ---
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

SatID = [1, 15, 30]
time_grid = np.linspace(0, 3, 1000)

fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
fig.suptitle("Broadcast vs Precise Orbit Residuals (RAC Frame)", fontsize=16, fontweight='bold')

rac_labels = ['Radial Error [m]', 'Along-Track Error [m]', 'Cross-Track Error [m]']

for prn in SatID:
    # Retrieve synchronized positions (mocked array alignment)
    pos_brdc = np.column_stack((np.cos(time_grid), np.sin(time_grid), np.ones_like(time_grid))) * 26560000.0
    
    # Inject typical broadcast ephemeris degradation (1-3 meters of systematic drift)
    error_injection = np.column_stack((np.sin(time_grid)*1.5, np.cos(time_grid)*2.0, np.sin(time_grid)*0.5))
    pos_sp3 = pos_brdc + error_injection
    
    err_radial = np.zeros(len(time_grid))
    err_along = np.zeros(len(time_grid))
    err_cross = np.zeros(len(time_grid))
    
    # Frame decomposition via instantaneous cross products
    for k in range(len(time_grid)):
        r_vec = pos_sp3[k, :]
        v_vec = np.array([-r_vec[1], r_vec[0], 0.0]) # Instantaneous orbital velocity approximation
        
        # Define local unit vectors
        u_rad = r_vec / np.linalg.norm(r_vec)
        u_cross = np.cross(r_vec, v_vec)
        u_cross /= np.linalg.norm(u_cross)
        u_along = np.cross(u_cross, u_rad)
        
        diff_cartesian = pos_brdc[k, :] - pos_sp3[k, :]
        
        err_radial[k] = np.dot(diff_cartesian, u_rad)
        err_along[k] = np.dot(diff_cartesian, u_along)
        err_cross[k] = np.dot(diff_cartesian, u_cross)
        
    # Plot components
    axes[0].plot(time_grid, err_radial, linewidth=1.5, label=f"PRN {prn}")
    axes[1].plot(time_grid, err_along, linewidth=1.5)
    axes[2].plot(time_grid, err_cross, linewidth=1.5)

for i, ax in enumerate(axes):
    ax.set_ylabel(rac_labels[i], fontsize=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    if i == 0:
        ax.legend(loc='upper right')

axes[2].set_xlabel("Time [hours]", fontsize=12)
plt.tight_layout()

save_path = os.path.join(output_dir, "orbit_accuracy_comparison_rac.png")
plt.savefig(save_path, dpi=300)
print(f"Comparison plot successfully saved to: {save_path}")
plt.show()