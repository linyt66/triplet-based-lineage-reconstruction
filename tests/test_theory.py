import math
import unittest

import _path  # noqa: F401

from triplet_lineage.theory import (
    admissible_triplet_error,
    binomial_success_probability,
    recursive_tmc_bound,
    sample_complexity_bound,
    single_level_oto_bound,
    solve_delta_star,
)


class TheoryTests(unittest.TestCase):
    def test_solve_delta_star_matches_closed_form(self):
        observed = solve_delta_star(d=0.5, lam=0.25, q=0.04)
        expected = 0.6 * (
            math.exp(-0.25 * 0.5) * (1 - 0.04)
            + 0.04 * math.exp(-0.25 * (2 - 0.5))
        )
        self.assertAlmostEqual(observed, expected)

    def test_sample_complexity_increases_with_missing_rate(self):
        delta_star = solve_delta_star(d=1.0, lam=0.25, q=0.04)
        no_missing = sample_complexity_bound(
            error=0.05,
            q=0.04,
            l_star=1 / 9,
            lam=0.25,
            delta_star=delta_star,
            p_miss=0.0,
        )
        with_missing = sample_complexity_bound(
            error=0.05,
            q=0.04,
            l_star=1 / 9,
            lam=0.25,
            delta_star=delta_star,
            p_miss=0.1,
        )
        self.assertGreater(with_missing, no_missing)

    def test_tmc_bound_decreases_with_triplet_error(self):
        self.assertGreater(recursive_tmc_bound(0.01), recursive_tmc_bound(0.1))
        self.assertGreater(single_level_oto_bound(0.01), single_level_oto_bound(0.1))

    def test_admissible_triplet_error_hits_target_bound(self):
        target = 0.6
        p_error = admissible_triplet_error(target, bound="recursive")

        self.assertGreaterEqual(recursive_tmc_bound(p_error), target - 1e-6)
        self.assertLess(recursive_tmc_bound(p_error + 1e-4), target)

    def test_binomial_success_probability_matches_simple_case(self):
        # With two triplets, target accuracy 1.0 requires both to be correct.
        self.assertAlmostEqual(
            binomial_success_probability(2, target_accuracy=1.0, p_error=0.2),
            0.64,
        )

    def test_invalid_parameters_raise_clear_errors(self):
        with self.assertRaises(ValueError):
            sample_complexity_bound(
                error=0.05,
                q=0.04,
                l_star=1 / 9,
                lam=0.25,
                delta_star=1.0,
                p_miss=1.0,
            )
        with self.assertRaises(ValueError):
            admissible_triplet_error(0.99)


if __name__ == "__main__":
    unittest.main()
