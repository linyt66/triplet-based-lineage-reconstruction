"""Theoretical bounds used in OTO/TMC parameter-scaling figures."""

from __future__ import annotations

import math

DEFAULT_MAXCUT_APPROXIMATION = 0.878


def solve_delta_star(d: float, lam: float, q: float) -> float:
    """Evaluate the depth-dependent separation term delta star."""
    _validate_nonnegative("d", d)
    _validate_positive("lam", lam)
    _validate_probability("q", q)
    return 0.6 * (
        math.exp(-lam * d) * (1 - q)
        + q * math.exp(-lam * (2 - d))
    )


def compute_delta_star(lam: float, q: float, d_max: float) -> float:
    """Compute optimal delta star by minimizing over d in [0, d_max]."""
    _validate_positive("lam", lam)
    _validate_probability("q", q)
    _validate_positive("d_max", d_max)

    from scipy.optimize import minimize_scalar

    result = minimize_scalar(
        solve_delta_star,
        bounds=(0, d_max),
        method="bounded",
        args=(lam, q),
    )
    return float(result.fun)


def sample_complexity_bound(
    error: float,
    q: float,
    l_star: float,
    lam: float,
    delta_star: float,
    p_miss: float = 0.0,
) -> float:
    """Compute the minimum site count implied by the manuscript bound."""
    _validate_probability_open("error", error)
    _validate_probability("q", q)
    _validate_positive("l_star", l_star)
    _validate_positive("lam", lam)
    _validate_positive("delta_star", delta_star)
    _validate_probability("p_miss", p_miss)
    if p_miss >= 1:
        raise ValueError("p_miss must be less than 1 for a finite bound.")

    log_term = -32 * math.log(error)

    k1 = log_term * q / (
        l_star**2 * delta_star**2 * (1 - p_miss) ** 5
    )
    k2 = log_term * (l_star + q * (1 - math.exp(-lam))) / (
        0.6
        * lam
        * l_star**2
        * delta_star
        * (1 - q + q * math.exp(-2 * lam))
        * (1 - p_miss) ** 3
    )

    return float(max(k1, k2))


def single_level_oto_bound(
    p_error: float,
    alpha_opt: float = DEFAULT_MAXCUT_APPROXIMATION,
) -> float:
    """Return the manuscript single-level OTO accuracy lower bound."""
    _validate_probability("p_error", p_error)
    _validate_probability_open("alpha_opt", alpha_opt)
    bound = (1 / 3 + (10 * alpha_opt - 6) / 9) - (2 * alpha_opt / 3) * p_error
    return _clip_unit_interval(bound)


def recursive_tmc_bound(
    p_error: float,
    alpha_opt: float = DEFAULT_MAXCUT_APPROXIMATION,
) -> float:
    """Return the manuscript recursive TMC accuracy lower bound."""
    _validate_probability("p_error", p_error)
    _validate_probability_open("alpha_opt", alpha_opt)
    bound = ((5 * alpha_opt / 3 - 2 / 3) - alpha_opt * p_error) * (1 - p_error)
    return _clip_unit_interval(bound)


def admissible_triplet_error(
    target_accuracy: float,
    alpha_opt: float = DEFAULT_MAXCUT_APPROXIMATION,
    bound: str = "recursive",
    tolerance: float = 1e-8,
) -> float:
    """Infer the largest triplet error allowed by a target tree accuracy."""
    _validate_probability("target_accuracy", target_accuracy)
    _validate_probability_open("alpha_opt", alpha_opt)
    _validate_positive("tolerance", tolerance)

    bound_fn = _select_bound(bound)
    if bound_fn(0.0, alpha_opt) < target_accuracy:
        raise ValueError("target_accuracy exceeds the selected bound at p_error=0.")

    lo = 0.0
    hi = 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if bound_fn(mid, alpha_opt) >= target_accuracy:
            lo = mid
        else:
            hi = mid
        if hi - lo <= tolerance:
            break
    return float(lo)


def binomial_success_probability(
    n_triplets: int,
    target_accuracy: float,
    p_error: float,
) -> float:
    """Approximate probability that empirical triplet accuracy hits a target."""
    if n_triplets < 1:
        raise ValueError("n_triplets must be positive.")
    _validate_probability("target_accuracy", target_accuracy)
    _validate_probability("p_error", p_error)

    min_successes = math.ceil(target_accuracy * n_triplets)
    p_correct = 1 - p_error
    try:
        from scipy.stats import binom

        return float(binom.sf(min_successes - 1, n_triplets, p_correct))
    except Exception:
        return _binomial_tail_exact(n_triplets, min_successes, p_correct)


def _binomial_tail_exact(n: int, k_min: int, p: float) -> float:
    if k_min <= 0:
        return 1.0
    if k_min > n:
        return 0.0
    return float(
        sum(
            math.comb(n, k) * p**k * (1 - p) ** (n - k)
            for k in range(k_min, n + 1)
        )
    )


def _select_bound(bound: str):
    if bound == "recursive":
        return recursive_tmc_bound
    if bound == "single_level":
        return single_level_oto_bound
    raise ValueError("bound must be 'recursive' or 'single_level'.")


def _validate_probability(name: str, value: float) -> None:
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must be between 0 and 1.")


def _validate_probability_open(name: str, value: float) -> None:
    if not 0 < value < 1:
        raise ValueError(f"{name} must be strictly between 0 and 1.")


def _validate_positive(name: str, value: float) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive.")


def _validate_nonnegative(name: str, value: float) -> None:
    if value < 0:
        raise ValueError(f"{name} must be nonnegative.")


def _clip_unit_interval(value: float) -> float:
    return float(min(1.0, max(0.0, value)))
