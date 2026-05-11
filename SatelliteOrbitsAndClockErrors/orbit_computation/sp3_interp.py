import numpy as np


def _lagrange_at(x_nodes, y_nodes, x):
    """Numerically stable direct Lagrange interpolation.

    x_nodes (N,), y_nodes (N,), x scalar.
    """
    N = len(x_nodes)
    total = 0.0
    for i in range(N):
        term = y_nodes[i]
        for j in range(N):
            if j != i:
                # Step-by-step division prevents intermediate float overflow
                term *= (x - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
        total += term
    return total


def interp_precise_orbits(
    pos, time, interp_time, order=10, convert_to_meters=True
):
    """Interpolate SP3 satellite positions to a dense time grid.

    Parameters
    ----------
    pos               : (N, 3) ECEF positions from SP3 [km]
    time              : (N,)   SP3 epochs [seconds]
    interp_time       : (M,)   target times [same units as `time`]
    order             : polynomial order (default 10 -> 11-point window)
    convert_to_meters : bool, multiplies output by 1000 to match broadcast units

    Returns
    -------
    orbit : (M, 3) interpolated positions [km or m]
    """
    pos = np.asarray(pos, dtype=float)
    time = np.asarray(time, dtype=float)
    target = np.atleast_1d(np.asarray(interp_time, dtype=float))

    # 1. STRICT GUARDRAIL: Prevent silent polynomial extrapolation
    if np.any(target < time[0]) or np.any(target > time[-1]):
        raise ValueError(
            f"Extrapolation detected! Target times (min: {target.min():.2f}, max: {target.max():.2f}) "
            f"fall outside the SP3 time range (min: {time[0]:.2f}, max: {time[-1]:.2f}).\n"
            "Verify that both time arrays use the exact same epoch origin and units."
        )

    n_nodes = order + 1
    orbit = np.empty((target.size, 3))

    for k, t_k in enumerate(target):
        # Centered window of n_nodes SP3 points around t_k
        i_center = int(np.searchsorted(time, t_k))
        i_start = max(0, i_center - n_nodes // 2)
        i_end = min(len(time), i_start + n_nodes)
        i_start = max(0, i_end - n_nodes)
        sl = slice(i_start, i_end)

        for d in range(3):
            orbit[k, d] = _lagrange_at(time[sl], pos[sl, d], t_k)

    # 2. UNIT ALIGNMENT: Convert km to m to match broadcast magnitude
    if convert_to_meters:
        orbit *= 1000.0

    return orbit