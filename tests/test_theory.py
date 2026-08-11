import math
import unittest

import _path  # noqa: F401

from triplet_lineage.theory import sample_complexity_bound, solve_delta_star


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


if __name__ == "__main__":
    unittest.main()
