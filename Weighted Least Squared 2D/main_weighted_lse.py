import numpy as np
import matplotlib.pyplot as plt
from error_ellipse import error_ellipse

# --- Integrated Data Generation ---
def data_generation(problem_type, example_number):
    x0_true = 4.4
    y0_true = 3.0
    ex = -2.0
    ey = 1.5
    R_pos = np.array([x0_true + ex, y0_true + ey], dtype=float)
    
    if example_number == 1:
        GS_x = np.array([9.0, 0.2, 3.7, -9.2, -5.5, 8.4, -6.0, 10.4, 10.5])
        GS_y = np.array([4.0, 7.7, 1.5, 5.8, -2.8, -5.7, 1.0, -0.8, 8.0])
        p_meas = np.array([4.7, 6.3, 1.7, 13.9, 11.5, 9.6, 10.6, 7.2, 7.9])
    elif example_number == 2:
        GS_x = np.array([9.8, 6.2, 6.7, 0.6, 3.3, -2.0, -2.8])
        GS_y = np.array([7.8, 6.7, 3.2, 2.5, -2.7, -1.0, -6.0])
        p_meas = np.array([7.2, 4.1, 2.3, 3.8, 5.8, 7.6, 11.5])
    else:
        raise ValueError("example_number must be 1 or 2")
        
    GS = np.column_stack((GS_x, GS_y))
    GS_n = len(GS_x)
    
    # Universally initialize states to prevent undefined variables
    ctau_true = 0.0
    R_ctau = 0.0
    Dx = np.zeros(3 if problem_type == 2 else 2)
    
    return x0_true, y0_true, GS, GS_n, p_meas, Dx, R_pos, ctau_true, R_ctau

# --- Configuration ---
problem_type = 2     # 1 = without clock bias; 2 = with receiver clock bias
example_number = 1   # 1 = good geometry; 2 = aligned stations
weighted = 1         # 0 = no weights; 1 = identity variance; 2 = diagonal approx

# Load Data
x0_true, y0_true, GS, GS_n, p_meas, Dx, R_pos, ctau_true, R_ctau = data_generation(problem_type, example_number)

# --- Initialization ---
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
    H = np.zeros((GS_n, 2))
    Dx_tot = np.zeros((k_max, 2))
    Cxx_tot = np.zeros((2, 2, k_max))
    Q_tot = np.zeros((2, 2, k_max))
elif problem_type == 2:
    H = np.zeros((GS_n, 3))
    estimated_ctau = np.zeros(k_max + 1)
    estimated_ctau[0] = R_ctau
    Dx_tot = np.zeros((k_max, 3))
    Cxx_tot = np.zeros((3, 3, k_max))
    Q_tot = np.zeros((3, 3, k_max))

# --- Main Iterative Loop ---
while stop_check > tol and k < k_max:
    # Distances and pseudoranges
    d_calc = np.sqrt((GS[:, 0] - R_pos[0])**2 + (GS[:, 1] - R_pos[1])**2)
    p_calc = d_calc + R_ctau
    
    # Pre-fit residuals
    Y_res = p_meas - p_calc
    r_norm_pre[k] = np.dot(Y_res, Y_res)
    mean_Y[k] = np.mean(Y_res)
    var_Y[k] = (1.0 / (GS_n - 1)) * np.dot(Y_res, Y_res)
    
    # Design matrix H formulation
    H[:, 0] = -(GS[:, 0] - R_pos[0]) / d_calc
    H[:, 1] = -(GS[:, 1] - R_pos[1]) / d_calc
    if problem_type == 2:
        H[:, 2] = 1.0
        
    # Weights and Covariance (Cpp) allocation
    if weighted == 0:
        W = np.eye(GS_n)
        Cpp = var_Y[k] * np.eye(GS_n)
    elif weighted == 1:
        Cpp = var_Y[k] * np.eye(GS_n)
        W = np.diag(1.0 / np.diagonal(Cpp))  # Fast inverse for diagonal matrices
    elif weighted == 2:
        Cpp = np.diag(Y_res**2)
        W = np.diag(1.0 / np.diagonal(Cpp))
        
    # Solve Weighted Normal Equations
    AA = H.T @ W @ H
    bb = H.T @ W @ Y_res
    Dx_new = np.linalg.solve(AA, bb)
    
    # Evaluate convergence on planimetric components
    stop_check = np.linalg.norm(Dx_new[0:2])
    
    # Update state vectors
    Dx = Dx_new
    R_pos = R_pos + Dx[0:2]
    Positions[k + 1, :] = R_pos
    Dx_tot[k, :] = Dx
    
    # Update clock bias
    if problem_type == 2:
        R_ctau = R_ctau + Dx[2]
        estimated_ctau[k + 1] = R_ctau
        
    # Evaluate post-fit residuals
    res_post = Y_res - H @ Dx
    r_norm_post[k] = np.dot(res_post, res_post)
    
    # Compute Cofactor (Q) and Covariance (Cxx) matrices
    Q = np.linalg.inv(H.T @ H)
    Q_tot[:, :, k] = Q
    
    if weighted == 0:
        Cxx = Q @ H.T @ Cpp @ H @ Q
    else:
        Cxx = np.linalg.inv(AA)
        
    Cxx_tot[:, :, k] = Cxx
    k += 1

# Truncate arrays to completed iterations
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

# --- Visualization Dashboard ---
# Chart 1: Planimetric Convergence
plt.figure(figsize=(7, 7))
plt.plot(Positions[:, 0], Positions[:, 1], 'k+--', label='Trajectory')
plt.plot(Positions[0, 0], Positions[0, 1], 'mo--', linewidth=2, label='Initial Guess')
plt.plot(Positions[-1, 0], Positions[-1, 1], 'bx--', linewidth=2, label='Final Estimate')
plt.plot(GS[:, 0], GS[:, 1], 'r*', label='Ground Stations')
plt.plot(x0_true, y0_true, 'g+', markeredgewidth=2, markersize=10, label='True Position')
plt.xlabel('East [cm]')
plt.ylabel('North [cm]')
plt.title('Planimetric Position Convergence')
plt.axis('equal')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

# Chart 2: Clock Bias Convergence
if problem_type == 2:
    plt.figure(figsize=(8, 4))
    plt.plot(range(k + 1), estimated_ctau, 'r--x', label='Estimated ctau')
    plt.plot(range(k + 1), ctau_true * np.ones(k + 1), 'b-', label='True ctau')
    plt.xlabel('Iteration')
    plt.ylabel('Clock Bias [cm]')
    plt.title('Receiver Clock Bias Estimation')
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend()
    plt.tight_layout()

# Chart 3: Position Uncertainties (Confidence Ellipses)
plt.figure(figsize=(7, 7))
plt.plot(Positions[:, 0], Positions[:, 1], 'k+--', label='Trajectory')
plt.plot(Positions[0, 0], Positions[0, 1], 'mo--', linewidth=2, label='Initial Guess')
plt.plot(Positions[-1, 0], Positions[-1, 1], 'bx--', linewidth=2, label='Final Estimate')
plt.plot(x0_true, y0_true, 'g+', markeredgewidth=2, markersize=10, label='True Position')

LevConf = 0.9889  # Mapping to k=3 (3-sigma)
for i in range(k):
    error_ellipse(Cxx_tot[0:2, 0:2, i], Positions[i + 1, :], conf=LevConf, color='dodgerblue')
    
plt.xlabel('East [cm]')
plt.ylabel('North [cm]')
plt.title('Positions with 98.89% Confidence Ellipses')
plt.axis('equal')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

# Chart 4: Residual Comparison
plt.figure(figsize=(8, 4))
plt.plot(range(k), r_norm_pre, 'r--o', label='Pre-fit Residual Norm')
plt.plot(range(1, k + 1), r_norm_post, 'b--x', label='Post-fit Residual Norm')
plt.title('Residual Norm Minimization Over Iterations')
plt.xlabel('Iteration')
plt.ylabel('Sum of Squared Residuals [cm²]')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

plt.show()