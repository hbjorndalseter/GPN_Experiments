import os
import numpy as np
import matplotlib.pyplot as plt

# --- Output Directory Setup ---
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# --- Reference Setup ---
# Ground Station Coordinates (Padova LLA: Lat, Lon in degrees, Alt in meters)
RecPos_lla = np.array([45.4108534, 11.8917094, 70.0413117])
maskAngle = 15.0  # Elevation Mask [degrees]
SatID = list(range(1, 13))  # PRN 1 through 12
time_grid_hours = np.linspace(0, 24, 500)

def ecef_to_enu_az_el(sat_ecef, rec_lat_deg, rec_lon_deg):
    """Transforms ECEF vector to local Azimuth and Elevation angles."""
    lat = np.deg2rad(rec_lat_deg)
    lon = np.deg2rad(rec_lon_deg)
    
    # Transformation Matrix: ECEF to Local ENU Frame
    R = np.array([
        [-np.sin(lon),              np.cos(lon),             0.0],
        [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)],
        [ np.cos(lat)*np.cos(lon),  np.cos(lat)*np.sin(lon), np.sin(lat)]
    ])
    
    enu = R @ sat_ecef
    e, n, u = enu[0], enu[1], enu[2]
    
    horizontal_dist = np.sqrt(e**2 + n**2)
    elevation = np.rad2deg(np.arctan2(u, horizontal_dist))
    azimuth = np.rad2deg(np.arctan2(e, n)) % 360.0
    
    return azimuth, elevation

# Pre-allocate visibility arrays
Az_matrix = np.full((len(time_grid_hours), len(SatID)), np.nan)
El_matrix = np.full((len(time_grid_hours), len(SatID)), np.nan)
Vis_matrix = np.zeros((len(time_grid_hours), len(SatID)), dtype=bool)

# Process visibility matrix over time
for idx, prn in enumerate(SatID):
    # Simulate orbital passes over the station
    elev_pass = 50.0 * np.sin(time_grid_hours * 0.5 + prn) + 10.0
    azim_pass = (time_grid_hours * 15.0 + prn * 30.0) % 360.0
    
    # Apply Masking Angle condition
    valid_mask = elev_pass >= maskAngle
    
    Az_matrix[valid_mask, idx] = azim_pass[valid_mask]
    El_matrix[valid_mask, idx] = elev_pass[valid_mask]
    Vis_matrix[:, idx] = valid_mask

# --- Dual Visualization Dashboard ---
fig = plt.figure(figsize=(14, 6))

# Subplot 1: Polar Skyplot
ax_polar = fig.add_subplot(121, projection='polar')
ax_polar.set_theta_zero_location('N')
ax_polar.set_theta_direction(-1)
ax_polar.set_ylim(0, 90)
ax_polar.set_yticks([0, 30, 60, 90])
ax_polar.set_yticklabels(['90°', '60°', '30°', '0°'])  # Inverted radial mapping
ax_polar.set_title(f"User Skyplot (Elevation Mask: {maskAngle}°)", fontsize=14, fontweight='bold', pad=20)

for idx, prn in enumerate(SatID):
    az_rad = np.deg2rad(Az_matrix[:, idx])
    el_mapped = 90.0 - El_matrix[:, idx]  # Center is 90 deg zenith
    
    if np.any(~np.isnan(el_mapped)):
        p = ax_polar.plot(az_rad, el_mapped, linewidth=1.5, label=f"PRN {prn}")
        # Mark instantaneous terminal positions
        valid_indices = np.where(~np.isnan(el_mapped))[0]
        last_idx = valid_indices[-1]
        ax_polar.scatter(az_rad[last_idx], el_mapped[last_idx], color=p[0].get_color(), marker='o')

# Subplot 2: Satellite Visibility Over Time
ax_vis = fig.add_subplot(122)
nSatView = np.sum(Vis_matrix, axis=1)

ax_vis.plot(time_grid_hours, nSatView, color='navy', linewidth=2)
ax_vis.axhline(np.mean(nSatView), color='red', linestyle='--', linewidth=1, label="Mean Availability")
ax_vis.set_title("Instantaneous Constellation Visibility", fontsize=14, fontweight='bold')
ax_vis.set_xlabel("Time [hours]", fontsize=12)
ax_vis.set_ylabel("Number of Visible Satellites", fontsize=12)
ax_vis.grid(True, linestyle=':', alpha=0.6)
ax_vis.legend()

plt.tight_layout()

save_path = os.path.join(output_dir, "skyplot_and_visibility.png")
plt.savefig(save_path, dpi=300)
print(f"Skyplot dashboard successfully saved to: {save_path}")
plt.show()