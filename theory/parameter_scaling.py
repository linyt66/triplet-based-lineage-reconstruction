import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfcinv
from matplotlib.cm import get_cmap
from plot_style import configure_journal_style
from matplotlib.ticker import ScalarFormatter

# =========================================================
# delta* function
# =========================================================
def solve_delta_star(d, lam, q):
    return 0.6 * (
        np.exp(-lam * d) * (1 - q)
        + q * np.exp(-lam * (2 - d))
    )

def compute_delta_star(lam, q, d_max):
    """Compute optimal delta* by minimizing over d ∈ [0, d_max]."""
    res = minimize_scalar(
        solve_delta_star,
        bounds=(0, d_max),
        method="bounded",
        args=(lam, q)
    )
    return res.fun


def sample_complexity_bound(
    error, q, l_star, lam, delta_star, p_miss=0.0
):
    """
    Compute minimal required number of sites n_site^min
    as max(k1, k2).
    """
    log_term = -32 * np.log(error)

    k1 = log_term * q / (l_star**2 * delta_star**2 * (1 - p_miss)**5)

    k2 = log_term * (l_star + q * (1 - np.exp(-lam))) / (
        0.6 * lam * l_star**2 * delta_star
        * (1 - q + q * np.exp(-2 * lam))
        * (1 - p_miss)**3
    )

    return max(k1, k2)

# Fixed parameters
# =========================================================
error = 0.05
q0 = 0.04
p_miss = 0.05
l_star0 = 1 / 9
d0 = 1
lambda0 = 0.25

# Precompute delta*
delta_star_0 = compute_delta_star(lambda0, q0, d0)

# =========================================================
# (c) error dependence
# =========================================================
error_values = np.logspace(np.log10(0.001), np.log10(0.3), 10)
k_values_error = np.array([
    sample_complexity_bound(
        err, q0, l_star0, lambda0, delta_star_0
    )
    for err in error_values
])

# =========================================================
# (d) l_min dependence
# =========================================================
l_star_values = np.linspace(1 / 20, 1 / 3, 13)
k_values_l_star = np.array([
    sample_complexity_bound(
        error, q0, l_star, lambda0, delta_star_0
    )
    for l_star in l_star_values
])

# =========================================================
# (e) lambda dependence
# =========================================================
lambda_values = np.logspace(-2, np.log10(5), 15)
k_values_lambda = []

for lam in lambda_values:
    delta_star = compute_delta_star(lam, q0, d0)
    k_values_lambda.append(
        sample_complexity_bound(
            error, q0, l_star0, lam, delta_star, p_miss
        )
    )
k_values_lambda = np.array(k_values_lambda)

# =========================================================
# Parameter grids
# =========================================================
q_values = np.logspace(np.log10(0.001), np.log10(0.8), 15)
d_values = np.linspace(0.1, 1, 15)
p_miss_values = np.logspace(np.log10(0.001), np.log10(0.5), 15)

# =========================================================
# k(q)
# =========================================================
k_values_q = np.array([
    sample_complexity_bound(
        error,
        q,
        l_star0,
        lambda0,
        compute_delta_star(lambda0, q, d0)
    )
    for q in q_values
])

# =========================================================
# k(d)
# =========================================================
delta_star_d = compute_delta_star(lambda0, q0, d0)
k_values_d = np.array([
    sample_complexity_bound(
        error,
        q0,
        l_star0,
        lambda0,
        compute_delta_star(lambda0, q0, d)
    )
    for d in d_values
])

# =========================================================
# k(p_miss)
# =========================================================
delta_star_0 = compute_delta_star(lambda0, q0, d0)
k_values_p_miss = np.array([
    sample_complexity_bound(
        error,
        q0,
        l_star0,
        lambda0,
        delta_star_0,
        p_miss
    )
    for p_miss in p_miss_values
])
