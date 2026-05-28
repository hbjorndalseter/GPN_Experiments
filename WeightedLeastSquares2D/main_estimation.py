import os
import numpy as np
import matplotlib.pyplot as plt
import data_generation

# --- Ensure Output Directory Exists ---
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# --- Configuration ---
problem_type = 1    # 1 = without the receiver clock bias; 2 = with clock bias
example_number = 1  # 1 = good geometry; 2 = aligned stations

# Load data
x0_true, y0_true, GS, GS_n, p_meas, Dx, R_pos, ctau_true, R_ctau = \
    data_generation.data_generation(problem_type, example_number)

# --- Loop Parameters Initialization ---
tol = 0.1          # Convergence threshold for position update norm
stop_check = 1.0   # Evaluated distance from previous position
k_max = 15         # Maximum number of iterations
k = 0              # Iteration counter

# Pre-allocate arrays
Positions = np.zeros((k_max + 1, 2))
Positions[0, :] = R_pos
r_norm_pre = np.zeros(k_max)   # Sum of squared pre-fit residuals

if problem_type == 1:
    n_states = 2
    H = np.zeros((GS_n, 2))
    Dx_tot = np.zeros((k_max, 2))
    header = f"{'Iter':>4} | {'East (x)':>14} | {'North (y)':>14} | {'||Dx||':>12} | {'Sum Sq Res (Pre)':>16}"
elif problem_type == 2:
    n_states = 3
    H = np.zeros((GS_n, 3))
    Dx_tot = np.zeros((k_max, 3))
    ctau_history = np.zeros(k_max + 1)
    ctau_history[0] = R_ctau
    header = (f"{'Iter':>4} | {'East (x)':>14} | {'North (y)':>14} | {'c*tau':>12} | "
              f"{'||Dx||':>12} | {'Sum Sq Res (Pre)':>16}")

print(header)
print("-" * len(header))

# --- Main Iterative Least-Squares Loop ---
while stop_check > tol and k < k_max:
    # Vectorized calculation of geometric distances
    d_calc = np.sqrt((GS[:, 0] - R_pos[0])**2 + (GS[:, 1] - R_pos[1])**2)

    if problem_type == 1:
        p_calc = d_calc.copy()
    elif problem_type == 2:
        p_calc = d_calc + R_ctau  # pseudorange = geometric range + clock bias

    # Pre-fit residuals: y_res = p_meas - p_calc
    Y_res = p_meas - p_calc
    r_norm_pre[k] = np.dot(Y_res, Y_res)

    # Vectorized construction of the Geometry / Design Matrix H (Jacobian)
    H[:, 0] = -(GS[:, 0] - R_pos[0]) / d_calc
    H[:, 1] = -(GS[:, 1] - R_pos[1]) / d_calc
    if problem_type == 2:
        H[:, 2] = 1.0  # partial derivative of pseudorange w.r.t. c*tau

    # Formulate Normal Equations: (H'H) * Dx = H'Y_res
    AA = np.dot(H.T, H)
    bb = np.dot(H.T, Y_res)

    # Solve system efficiently without direct inversion
    Dx_new = np.linalg.solve(AA, bb)

    # Evaluate stop condition (L2 norm of the state update vector)
    stop_check = np.linalg.norm(Dx_new)

    # Update receiver state
    Dx = Dx_new
    R_pos = R_pos + Dx[:2]
    if problem_type == 2:
        R_ctau = R_ctau + Dx[2]

    # Store variables
    Positions[k + 1, :] = R_pos
    Dx_tot[k, :] = Dx
    if problem_type == 2:
        ctau_history[k + 1] = R_ctau

    # Print iteration info
    if problem_type == 1:
        print(f" {k+1:3d} |  {R_pos[0]:12.5f}  |  {R_pos[1]:12.5f}  |  "
              f"{stop_check:10.5f}  |  {r_norm_pre[k]:14.5f}")
    elif problem_type == 2:
        print(f" {k+1:3d} |  {R_pos[0]:12.5f}  |  {R_pos[1]:12.5f}  |  "
              f"{R_ctau:10.5f}  |  {stop_check:10.5f}  |  {r_norm_pre[k]:14.5f}")

    k += 1

# Slice pre-allocated arrays down to actual completed iterations
Positions = Positions[:k + 1, :]
Dx_tot = Dx_tot[:k, :]
r_norm_pre = r_norm_pre[:k]
if problem_type == 2:
    ctau_history = ctau_history[:k + 1]

# --- Final Summary ---
print("\n--- Final Results ---")
print(f"True position:      ({x0_true:.4f}, {y0_true:.4f})")
print(f"Estimated position: ({R_pos[0]:.4f}, {R_pos[1]:.4f})")
pos_error = np.sqrt((R_pos[0] - x0_true)**2 + (R_pos[1] - y0_true)**2)
print(f"Position error:     {pos_error:.6f}")
if problem_type == 2:
    print(f"True clock bias:    {ctau_true:.4f}")
    print(f"Estimated c*tau:    {R_ctau:.4f}")
    print(f"Clock bias error:   {abs(R_ctau - ctau_true):.6f}")
print(f"Converged in {k} iterations")

# --- Visualization: Planimetric Network ---
prob_label = "with clock bias" if problem_type == 2 else "without clock bias"
plt.figure(figsize=(8, 8))
plt.title(f"Planimetric Network Convergence\n(Example {example_number}, {prob_label})",
          fontsize=14, fontweight='bold')

# Plot line-of-sight vectors from estimated final position to stations
for i in range(GS_n):
    plt.plot([R_pos[0], GS[i, 0]], [R_pos[1], GS[i, 1]],
             color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    plt.text(GS[i, 0] + 0.3, GS[i, 1] + 0.3, f"{i+1}",
             fontsize=10, fontweight='bold', color='navy')

# Plot ground stations
plt.scatter(GS[:, 0], GS[:, 1], color='blue', s=60, label='Ground Stations', zorder=3)

# Plot iteration trajectory
plt.plot(Positions[:, 0], Positions[:, 1], color='orange', marker='o',
         linestyle='-', linewidth=2, label='Iteration Path', zorder=4)
plt.scatter(Positions[0, 0], Positions[0, 1], color='purple', s=80,
            marker='X', label='Initial Guess', zorder=5)

# Plot true vs final estimated position
plt.scatter(x0_true, y0_true, color='green', s=100, marker='*',
            label='True Position', zorder=5)
plt.scatter(R_pos[0], R_pos[1], color='red', s=80,
            label='Estimated Position', zorder=5)

plt.xlabel("East [cm]", fontsize=12)
plt.ylabel("North [cm]", fontsize=12)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(loc='best')
plt.axis('equal')
plt.tight_layout()

# Save plot
tag = "biased" if problem_type == 2 else "unweighted"
save_path = os.path.join(output_dir, f"{tag}_convergence_ex{example_number}.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"\nPlot saved to: {save_path}")
plt.show()