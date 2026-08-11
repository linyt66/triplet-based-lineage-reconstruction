import unittest

import _path  # noqa: F401


class PackageImportTests(unittest.TestCase):
    def test_top_level_import_is_lightweight(self):
        import triplet_lineage

        self.assertIn("solve_delta_star", triplet_lineage.__all__)
        self.assertEqual(triplet_lineage.__name__, "triplet_lineage")

    def test_unknown_lazy_export_raises_attribute_error(self):
        import triplet_lineage

        with self.assertRaises(AttributeError):
            getattr(triplet_lineage, "does_not_exist")


if __name__ == "__main__":
    unittest.main()
