import unittest

import _path  # noqa: F401

from triplet_lineage.plotting import COLORS, log10_formatter


class PlottingTests(unittest.TestCase):
    def test_log10_formatter(self):
        self.assertEqual(log10_formatter(100), "2")
        self.assertEqual(log10_formatter(1), "0")
        self.assertEqual(log10_formatter(0), "")

    def test_palette_contains_required_series_colors(self):
        for color_name in ("blue", "orange", "green", "red", "purple", "dark_gray"):
            self.assertIn(color_name, COLORS)
            self.assertRegex(COLORS[color_name], r"^#[0-9a-fA-F]{6}$")


if __name__ == "__main__":
    unittest.main()
