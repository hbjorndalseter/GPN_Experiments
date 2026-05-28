import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2


def error_ellipse(C, mu, conf=0.9889, ax=None, **kwargs):
    """
    Plots an error ellipse defining the confidence region for a 2x2 covariance matrix.

    Parameters:
        C (ndarray): 2x2 Covariance matrix.
        mu (ndarray): 1D array representing the center (x, y) of the ellipse.
        conf (float): Confidence interval (between 0 and 1). Default is 0.9889 (k=3).
        ax (matplotlib.axes.Axes): Axes to plot on. Defaults to current axes.
        **kwargs: Additional keyword arguments passed to plt.plot.
    """
    if ax is None:
        ax = plt.gca()

    C = np.asarray(C)
    mu = np.asarray(mu).flatten()

    if C.shape != (2, 2):
        raise ValueError("Covariance matrix C must be exactly 2x2.")

    eigvals, eigvecs = np.linalg.eigh(C)

    if np.any(eigvals <= 0):
        raise ValueError("The covariance matrix must be positive definite.")

    k = np.sqrt(chi2.ppf(conf, df=2))

    theta = np.linspace(0, 2 * np.pi, 100)
    circle = np.column_stack((np.cos(theta), np.sin(theta)))

    ellipse = circle @ np.diag(np.sqrt(eigvals)) @ eigvecs.T

    x = mu[0] + k * ellipse[:, 0]
    y = mu[1] + k * ellipse[:, 1]

    if not kwargs:
        kwargs = {'color': 'blue', 'linestyle': '-', 'linewidth': 1, 'alpha': 0.6}

    ax.plot(x, y, **kwargs)