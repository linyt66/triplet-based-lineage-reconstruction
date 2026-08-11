"""Theoretical bounds used in OTO/TMC parameter-scaling figures."""

import math


def solve_delta_star(d: float, lam: float, q: float) -> float:
    """Evaluate the depth-dependent separation term delta star."""
    return 0.6 * (
        math.exp(-lam * d) * (1 - q)
        + q * math.exp(-lam * (2 - d))
    )


def compute_delta_star(lam: float, q: float, d_max: float) -> float:
    """Compute optimal delta star by minimizing over d in [0, d_max]."""
    from scipy.optimize import minimize_scalar

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
    log_term = -32 * math.log(error)

    k1 = log_term * q / (
        l_star**2 * delta_star**2 * (1 - p_miss)**5
    )
    k2 = log_term * (l_star + q * (1 - math.exp(-lam))) / (
        0.6
        * lam
        * l_star**2
        * delta_star
        * (1 - q + q * math.exp(-2 * lam))
        * (1 - p_miss)**3
    )

    return max(k1, k2)
