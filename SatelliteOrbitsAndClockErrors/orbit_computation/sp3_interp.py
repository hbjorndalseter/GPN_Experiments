"""SP3 precise orbit interpolation via centered Lagrange polynomials.

IGS standard practice: 9th–11th order polynomial across ~10–12 SP3 samples
centered around the target time. Cubic splines (your previous code) are
fine visually but introduce ~10 cm errors at window edges — Lagrange of
order 10 gets you below 1 cm interior to the day.
"""
import numpy as np


def _lagrange_at(x_nodes, y_nodes, x):
    """Direct Lagrange interpolation. x_nodes (N,), y_nodes (N,), x scalar."""
    N = len(x_nodes)
    total = 0.0
    for i in range(N):
        num = den = 1.0
        for j in range(N):
            if j != i:
                num *= (x - x_nodes[j])
                den *= (x_nodes[i] - x_nodes[j])
        total += y_nodes[i] * num / den
    return total


def interp_precise_orbits(pos, time, interp_time, order=10):
    """Interpolate SP3 satellite positions to a dense time grid.

    Parameters
    ----------
    pos         : (N, 3) ECEF positions from SP3   [km]
    time        : (N,)   SP3 epochs               [seconds, any consistent epoch]
    interp_time : (M,)   target times             [same units as `time`]
    order       : polynomial order (default 10 -> 11-point window).

    Returns
    -------
    orbit : (M, 3) interpolated positions [km]
    """
    pos  = np.asarray(pos, dtype=float)
    time = np.asarray(time, dtype=float)
    target = np.atleast_1d(np.asarray(interp_time, dtype=float))
    n_nodes = order + 1

    orbit = np.empty((target.size, 3))
    for k, t_k in enumerate(target):
        # Centered window of n_nodes SP3 points around t_k, clipped to bounds
        i_center = int(np.searchsorted(time, t_k))
        i_start  = max(0, i_center - n_nodes // 2)
        i_end    = min(len(time), i_start + n_nodes)
        i_start  = max(0, i_end - n_nodes)
        sl = slice(i_start, i_end)
        for d in range(3):
            orbit[k, d] = _lagrange_at(time[sl], pos[sl, d], t_k)
    return orbit