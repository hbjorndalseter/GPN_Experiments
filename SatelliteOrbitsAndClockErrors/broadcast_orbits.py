import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import gnss_lib_py as glp

# Load data transmitted with the Broadcast Navigation Message
# Automatically standardizes timing to 'gps_millis' and labels constellations/PRNs
filename = 'BRDC00IGS_R_20240310000_01D_MN.rnx'
rinex_nav = glp.RinexNav(filename)

# Extract broadcast ephemeris parameters for a specific satellite (e.g., GPS PRN 1)
# This exposes underlying arrays (M_0, e, sqrtA, SVclockBias, etc.) for propagation
prn_1_nav = rinex_nav.where("gnss_id", "gps")
if prn_1_nav is not None:
    prn_1_nav = prn_1_nav.where("sv_id", 1)

# --- Output Directory Setup ---
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# --- Configuration ---
filename = 'BRDC00IGS_R_20240310000_01D_MN.rnx'
constID = "GPS"  # Can be "GPS" or "GLONASS"
SatID = [1, 15, 30]
vector_time = np.linspace(0, 3, 1000)  # Time grid in hours
toe_start = 259200.0  # Reference Time of Ephemeris (seconds of week)

# --- Physical Constants (WGS84 / PZ-90) ---
MU_EARTH = 3.986004418e14  # m^3/s^2
OMEGA_EARTH = 7.292115e-5  # rad/s
J2_EARTH = 1.0826257e-3
RADIUS_EARTH = 6378137.0   # m

def compute_keplerian_gps(t_sec, toe, a, e, i0, omega, omega_dot, w, m0, delta_n):
    """Solves Kepler's equation and computes ECEF coordinates for GPS/Galileo."""
    n0 = np.sqrt(MU_EARTH / (a**3))
    n = n0 + delta_n
    tk = t_sec - toe
    mk = m0 + n * tk
    
    # Iteratively solve Kepler's Equation: E - e*sin(E) = M
    ek = mk
    for _ in range(10):
        ek = mk + e * np.sin(ek)
        
    vk = np.arctan2(np.sqrt(1.0 - e**2) * np.sin(ek), np.cos(ek) - e)
    phik = vk + w
    
    # Unperturbed radius and inclination (simplified for baseline propagation)
    rk = a * (1.0 - e * np.cos(ek))
    ik = i0
    
    xk_prime = rk * np.cos(phik)
    yk_prime = rk * np.sin(phik)
    
    omegak = omega + (omega_dot - OMEGA_EARTH) * tk - OMEGA_EARTH * toe
    
    x = xk_prime * np.cos(omegak) - yk_prime * np.cos(ik) * np.sin(omegak)
    y = xk_prime * np.sin(omegak) + yk_prime * np.cos(ik) * np.cos(omegak)
    z = yk_prime * np.sin(ik)
    
    return np.array([x, y, z]) / 1000.0  # Return in km

def glonass_derivatives(t, y_state):
    """Coupled GLONASS ECEF differential equations of motion."""
    r = np.sqrt(y_state[0]**2 + y_state[2]**2 + y_state[4]**2)
    r3 = r**3
    r5 = r**5
    
    # Common J2 perturbation coefficient
    j2_term = (3.0 / 2.0) * J2_EARTH * MU_EARTH * (RADIUS_EARTH**2) / r5
    z_term = 1.0 - (5.0 * y_state[4]**2) / (r**2)
    
    # Equations matching course workbook specifications
    dx = y_state[1]
    ddx = - (MU_EARTH * y_state[0]) / r3 + j2_term * y_state[0] * z_term + (OMEGA_EARTH**2) * y_state[0] + 2.0 * OMEGA_EARTH * y_state[3]
    
    dy = y_state[3]
    ddy = - (MU_EARTH * y_state[2]) / r3 + j2_term * y_state[2] * z_term + (OMEGA_EARTH**2) * y_state[2] - 2.0 * OMEGA_EARTH * y_state[1]
    
    dz = y_state[5]
    ddz = - (MU_EARTH * y_state[4]) / r3 + j2_term * y_state[4] * (3.0 - (5.0 * y_state[4]**2) / (r**2))
    
    return [dx, ddx, dy, ddy, dz, ddz]

# --- Main Orbit Computation Structure ---
SAT = {}
t_seconds = toe_start + vector_time * 3600.0

for prn in SatID:
    orbit = np.zeros((len(vector_time), 3))
    
    if constID == "GPS":
        # Mock broadcast parameters for demonstration (Semi-major axis ~26560 km)
        a_geo = 26560000.0
        ecc = 0.005
        inc = np.deg2rad(55.0)
        
        for j, t_val in enumerate(t_seconds):
            orbit[j, :] = compute_keplerian_gps(t_val, toe_start, a_geo, ecc, inc, 0.0, 0.0, 0.0, 0.0, 0.0)
            
    elif constID == "GLONASS":
        # Initial GLONASS state vector [x, vx, y, vy, z, vz] in meters and m/s
        y0 = [1.9e7, -2000.0, 1.5e7, 3000.0, 1.0e7, 1500.0]
        t_span = (0.0, vector_time[-1] * 3600.0)
        t_eval = vector_time * 3600.0
        
        sol = solve_ivp(glonass_derivatives, t_span, y0, t_eval=t_eval, method='RK45', rtol=1e-9, atol=1e-12)
        orbit = sol.y[0:6:2, :].T / 1000.0  # Extract x, y, z positions in km
        
    SAT[f"PRN_{prn}_brdc"] = orbit

# --- Visualization Dashboard ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_title(f"Broadcast Orbits Propagation ({constID})", fontsize=14, fontweight='bold')

# Render Earth Reference Sphere
u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
x_e = (RADIUS_EARTH / 1000.0) * np.cos(u) * np.sin(v)
y_e = (RADIUS_EARTH / 1000.0) * np.sin(u) * np.sin(v)
z_e = (RADIUS_EARTH / 1000.0) * np.cos(v)
ax.plot_wireframe(x_e, y_e, z_e, color='lightblue', alpha=0.3)

for prn in SatID:
    traj = SAT[f"PRN_{prn}_brdc"]
    ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], linewidth=2, label=f"Sat {prn}")
    ax.scatter(traj[0, 0], traj[0, 1], traj[0, 2], color='red', marker='o') # Start position

ax.set_xlabel("ECEF X [km]")
ax.set_ylabel("ECEF Y [km]")
ax.set_zlabel("ECEF Z [km]")
ax.legend()
plt.tight_layout()

save_path = os.path.join(output_dir, f"broadcast_propagation_{constID}.png")
plt.savefig(save_path, dpi=300)
print(f"Broadcast orbits successfully saved to: {save_path}")
plt.show()