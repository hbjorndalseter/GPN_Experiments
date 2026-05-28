import numpy as np


def data_generation(problem_type, example_number):
    """
    Generates the true position, initial guess, ground station coordinates,
    and measured pseudoranges for the 2D GNSS positioning problem.

    problem_type = 1: position only (x, y)
    problem_type = 2: position + receiver clock bias (x, y, c*tau)
    """
    # Receiver true position
    x0_true = 4.4
    y0_true = 3.0

    # Initial position error
    ex = -2.0
    ey = 1.5

    # Initial receiver position
    R_x0 = x0_true + ex
    R_y0 = y0_true + ey
    R_pos = np.array([R_x0, R_y0], dtype=float)

    # Select Ground Station layout based on example_number
    if example_number == 1:
        # Good positioning geometry (well-distributed around the receiver)
        GS_x = np.array([9.0, 0.2, 3.7, -9.2, -5.5, 8.4, -6.0, 10.4, 10.5])
        GS_y = np.array([4.0, 7.7, 1.5, 5.8, -2.8, -5.7, 1.0, -0.8, 8.0])
        p_meas = np.array([4.7, 6.3, 1.7, 13.9, 11.5, 9.6, 10.6, 7.2, 7.9])
    elif example_number == 2:
        # Poor positioning geometry (stations aligned)
        GS_x = np.array([9.8, 6.2, 6.7, 0.6, 3.3, -2.0, -2.8])
        GS_y = np.array([7.8, 6.7, 3.2, 2.5, -2.7, -1.0, -6.0])
        p_meas = np.array([7.2, 4.1, 2.3, 3.8, 5.8, 7.6, 11.5])
    else:
        raise ValueError("example_number must be 1 or 2")

    # Matrix of ground station positions
    GS = np.column_stack((GS_x, GS_y))
    GS_n = len(GS_x)

    # Initialize biases and state variables based on problem type
    if problem_type == 1:
        # No clock bias
        Dx = np.zeros(2)
        ctau_true = 0.0
        R_ctau = 0.0
    elif problem_type == 2:
        # With receiver clock bias
        Dx = np.zeros(3)          # state vector: [dx, dy, d(c*tau)]
        ctau_true = 2.0           # true clock bias (same units as pseudorange)
        R_ctau = 0.0              # initial guess for clock bias
        # Add the clock bias to the measured pseudoranges
        p_meas = p_meas + ctau_true
    else:
        raise ValueError("problem_type must be 1 or 2")

    return x0_true, y0_true, GS, GS_n, p_meas, Dx, R_pos, ctau_true, R_ctau