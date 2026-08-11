import importlib.util
import unittest

import _path  # noqa: F401


@unittest.skipUnless(
    importlib.util.find_spec("networkx")
    and importlib.util.find_spec("numpy")
    and importlib.util.find_spec("pandas"),
    "networkx, numpy, and pandas are required for reconstruction graph tests",
)
class ReconstructionGraphTests(unittest.TestCase):
    def test_construct_triplet_graph_accumulates_weights(self):
        from triplet_lineage.reconstruction import construct_triplet_graph

        graph = construct_triplet_graph([
            ("a", "b", "c"),
            ("a", "b", "d"),
        ])

        self.assertEqual(graph["a"]["b"]["weight"], -4)
        self.assertEqual(graph["a"]["c"]["weight"], 1)
        self.assertEqual(graph["b"]["c"]["weight"], 1)
        self.assertEqual(graph["a"]["d"]["weight"], 1)
        self.assertEqual(graph["b"]["d"]["weight"], 1)

    def test_infer_triplets_from_minimum_hamming_distance(self):
        import pandas as pd
        from triplet_lineage.reconstruction import infer_triplets_from_mutations

        class FakeCassiopeiaTree:
            leaves = ["a", "b", "c"]
            character_matrix = pd.DataFrame(
                [[1, 0, 0], [1, 0, 1], [0, 1, 1]],
                index=leaves,
            )

        self.assertEqual(
            infer_triplets_from_mutations(FakeCassiopeiaTree()),
            [("a", "b", "c")],
        )


if __name__ == "__main__":
    unittest.main()
