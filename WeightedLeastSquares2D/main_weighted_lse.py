import os
import numpy as np
import matplotlib.pyplot as plt
from error_ellipse import error_ellipse
import data_generation

# --- Ensure Output Directory Exists ---
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# --- Configuration ---
problem_type = 2     # 1 = without clock bias; 2 = with receiver clock bias
example_number = 2   # 1 = good geometry; 2 = aligned stations
weighted = 2       # 0 = no weights; 1 = identity variance; 2 = diagonal approx

# Load data from shared module
x0_true, y0_true, GS, GS_n, p_meas, Dx, R_pos, ctau_true, R_ctau = \
    data_generation.data_generation(problem_type, example_number)

# --- Loop Parameters ---
tol = 0.1
stop_check = 1.0
k_max = 15
k = 0

Positions = np.zeros((k_max + 1, 2))
Positions[0, :] = R_pos

r_norm_pre = np.zeros(k_max)
r_norm_post = np.zeros(k_max)
mean_Y = np.zeros(k_max)
var_Y = np.zeros(k_max)

if problem_type == 1:
    n_states = 2
    H = np.zeros((GS_n, 2))
    Dx_tot = np.zeros((k_max, 2))
    Cxx_tot = np.zeros((2, 2, k_max))
    Q_tot = np.zeros((2, 2, k_max))
elif problem_type == 2:
    n_states = 3
    H = np.zeros((GS_n, 3))
    estimated_ctau = np.zeros(k_max + 1)
    estimated_ctau[0] = R_ctau
    Dx_tot = np.zeros((k_max, 3))
    Cxx_tot = np.zeros((3, 3, k_max))
    Q_tot = np.zeros((3, 3, k_max))

# --- Header ---
prob_label = "with clock bias" if problem_type == 2 else "no clock bias"
weight_labels = {0: "unweighted", 1: "identity variance", 2: "diagonal approx"}
print(f"WLSE — Problem type {problem_type} ({prob_label}), "
      f"Example {example_number}, Weighting: {weight_labels[weighted]}")
print("=" * 90)

if problem_type == 1:
    header = (f"{'Iter':>4} | {'East (x)':>12} | {'North (y)':>12} | "
              f"{'||Dx||':>10} | {'PreFit SSR':>12} | {'PostFit SSR':>12}")
else:
    header = (f"{'Iter':>4} | {'East (x)':>12} | {'North (y)':>12} | {'c*tau':>10} | "
              f"{'||Dx||':>10} | {'PreFit SSR':>12} | {'PostFit SSR':>12}")
print(header)
print("-" * len(header))

# --- Main Iterative Weighted Least-Squares Loop ---
while stop_check > tol and k < k_max:
    # Geometric distances
    d_calc = np.sqrt((GS[:, 0] - R_pos[0])**2 + (GS[:, 1] - R_pos[1])**2)

    # Predicted pseudoranges
    if problem_type == 1:
        p_calc = d_calc.copy()
    else:
        p_calc = d_calc + R_ctau

    # Pre-fit residuals
    Y_res = p_meas - p_calc
    r_norm_pre[k] = np.dot(Y_res, Y_res)
    mean_Y[k] = np.mean(Y_res)
    var_Y[k] = (1.0 / (GS_n - 1)) * np.dot(Y_res, Y_res)

    # Design matrix (Jacobian)
    H[:, 0] = -(GS[:, 0] - R_pos[0]) / d_calc
    H[:, 1] = -(GS[:, 1] - R_pos[1]) / d_calc
    if problem_type == 2:
        H[:, 2] = 1.0

    # Construct weight matrix
    if weighted == 0:
        W = np.eye(GS_n)
        Cpp = var_Y[k] * np.eye(GS_n)
    elif weighted == 1:
        Cpp = var_Y[k] * np.eye(GS_n)
        W = np.diag(1.0 / np.diagonal(Cpp))
    elif weighted == 2:
        Cpp = np.diag(Y_res**2)
        W = np.diag(1.0 / np.diagonal(Cpp))

    # Weighted normal equations
    AA = H.T @ W @ H
    bb = H.T @ W @ Y_res
    Dx_new = np.linalg.solve(AA, bb)

    # Convergence check (position components only)
    stop_check = np.linalg.norm(Dx_new[0:2])

    # Update state
    Dx = Dx_new
    R_pos = R_pos + Dx[0:2]
    Positions[k + 1, :] = R_pos
    Dx_tot[k, :] = Dx

    if problem_type == 2:
        R_ctau = R_ctau + Dx[2]
        estimated_ctau[k + 1] = R_ctau

    # Post-fit residuals
    res_post = Y_res - H @ Dx
    r_norm_post[k] = np.dot(res_post, res_post)

    # Cofactor matrix (unweighted geometry)
    Q = np.linalg.inv(H.T @ H)
    Q_tot[:, :, k] = Q

    # Covariance of the estimated parameters
    if weighted == 0:
        Cxx = Q @ H.T @ Cpp @ H @ Q
    else:
        Cxx = np.linalg.inv(AA)

    Cxx_tot[:, :, k] = Cxx

    # Print iteration
    if problem_type == 1:
        print(f" {k+1:3d} |  {R_pos[0]:10.5f}  |  {R_pos[1]:10.5f}  |  "
              f"{stop_check:8.5f}  |  {r_norm_pre[k]:10.5f}  |  {r_norm_post[k]:10.5f}")
    else:
        print(f" {k+1:3d} |  {R_pos[0]:10.5f}  |  {R_pos[1]:10.5f}  |  {R_ctau:8.5f}  |  "
              f"{stop_check:8.5f}  |  {r_norm_pre[k]:10.5f}  |  {r_norm_post[k]:10.5f}")

    k += 1

# --- Trim arrays ---
Positions = Positions[0:k + 1, :]
r_norm_pre = r_norm_pre[0:k]
r_norm_post = r_norm_post[0:k]
mean_Y = mean_Y[0:k]
var_Y = var_Y[0:k]
Dx_tot = Dx_tot[0:k, :]
Q_tot = Q_tot[:, :, 0:k]
Cxx_tot = Cxx_tot[:, :, 0:k]
if problem_type == 2:
    estimated_ctau = estimated_ctau[0:k + 1]

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
print(f"\nFinal Cxx (position block):")
print(Cxx_tot[0:2, 0:2, -1])

# ====================================================================
# VISUALIZATION
# ====================================================================
tag = f"wlse_p{problem_type}_ex{example_number}_w{weighted}"

# Chart 1: Planimetric Convergence
fig1 = plt.figure(figsize=(7, 7))
plt.plot(Positions[:, 0], Positions[:, 1], 'k+--', label='Trajectory')
plt.plot(Positions[0, 0], Positions[0, 1], 'mo--', linewidth=2, label='Initial Guess')
plt.plot(Positions[-1, 0], Positions[-1, 1], 'bx--', linewidth=2, label='Final Estimate')
plt.plot(GS[:, 0], GS[:, 1], 'r*', markersize=10, label='Ground Stations')
plt.plot(x0_true, y0_true, 'g+', markeredgewidth=2, markersize=12, label='True Position')
for i in range(GS_n):
    plt.text(GS[i, 0] + 0.3, GS[i, 1] + 0.3, f"{i+1}",
             fontsize=9, fontweight='bold', color='navy')
plt.xlabel('East [cm]')
plt.ylabel('North [cm]')
plt.title(f'Planimetric Convergence — {prob_label}, Ex {example_number}')
plt.axis('equal')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()
fig1.savefig(os.path.join(output_dir, f"{tag}_planimetric.png"), dpi=300)

# Chart 2: Clock Bias Convergence (problem_type == 2 only)
if problem_type == 2:
    fig2 = plt.figure(figsize=(8, 4))
    plt.plot(range(k + 1), estimated_ctau, 'r--x', label='Estimated c·τ')
    plt.axhline(ctau_true, color='b', linestyle='-', label=f'True c·τ = {ctau_true}')
    plt.xlabel('Iteration')
    plt.ylabel('Clock Bias [cm]')
    plt.title(f'Receiver Clock Bias Estimation — Ex {example_number}')
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    fig2.savefig(os.path.join(output_dir, f"{tag}_clockbias.png"), dpi=300)

# Chart 3: Confidence Ellipses
fig3 = plt.figure(figsize=(7, 7))
plt.plot(Positions[:, 0], Positions[:, 1], 'k+--', label='Trajectory')
plt.plot(Positions[0, 0], Positions[0, 1], 'mo--', linewidth=2, label='Initial Guess')
plt.plot(Positions[-1, 0], Positions[-1, 1], 'bx--', linewidth=2, label='Final Estimate')
plt.plot(x0_true, y0_true, 'g+', markeredgewidth=2, markersize=12, label='True Position')

LevConf = 0.9889  # 3-sigma for 2 DOF
for i in range(k):
    error_ellipse(Cxx_tot[0:2, 0:2, i], Positions[i + 1, :],
                  conf=LevConf, color='dodgerblue')

plt.xlabel('East [cm]')
plt.ylabel('North [cm]')
plt.title(f'Positions with 98.89% Confidence Ellipses — Ex {example_number}')
plt.axis('equal')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()
fig3.savefig(os.path.join(output_dir, f"{tag}_ellipses.png"), dpi=300)

# Chart 4: Residual Comparison
fig4 = plt.figure(figsize=(8, 4))
iters = np.arange(1, k + 1)
plt.plot(iters, r_norm_pre, 'r--o', label='Pre-fit SSR')
plt.plot(iters, r_norm_post, 'b--x', label='Post-fit SSR')
plt.title(f'Residual Norm Minimization — Ex {example_number}')
plt.xlabel('Iteration')
plt.ylabel('Sum of Squared Residuals [cm²]')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()
fig4.savefig(os.path.join(output_dir, f"{tag}_residuals.png"), dpi=300)

print(f"\nAll plots saved to '{output_dir}/'")
plt.show()