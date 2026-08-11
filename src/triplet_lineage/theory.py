"""Theoretical bounds used in OTO/TMC parameter-scaling figures."""

import numpy as np
from scipy.optimize import minimize_scalar


def solve_delta_star(d: float, lam: float, q: float) -> float:
    """Evaluate the depth-dependent separation term delta star."""
    return 0.6 * (
        np.exp(-lam * d) * (1 - q)
        + q * np.exp(-lam * (2 - d))
    )


def compute_delta_star(lam: float, q: float, d_max: float) -> float:
    """Compute optimal delta star by minimizing over d in [0, d_max]."""
    result = minimize_scalar(
        solve_delta_star,
        bounds=(0, d_max),
        method="bounded",
        args=(lam, q),
    )
    return result.fun


def sample_complexity_bound(
    error: float,
    q: float,
    l_star: float,
    lam: float,
    delta_star: float,
    p_miss: float = 0.0,
) -> float:
    """Compute the minimum site count implied by the manuscript bound."""
    log_term = -32 * np.log(error)

    k1 = log_term * q / (
        l_star**2 * delta_star**2 * (1 - p_miss)**5
    )
    k2 = log_term * (l_star + q * (1 - np.exp(-lam))) / (
        0.6
        * lam
        * l_star**2
        * delta_star
        * (1 - q + q * np.exp(-2 * lam))
        * (1 - p_miss)**3
    )

    return max(k1, k2)
