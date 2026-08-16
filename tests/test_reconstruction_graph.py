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

    def test_construct_triplet_graph_accepts_confidence_weights(self):
        from triplet_lineage.reconstruction import construct_triplet_graph

        graph = construct_triplet_graph([("a", "b", "c", 0.5)])

        self.assertEqual(graph["a"]["b"]["weight"], -1)
        self.assertEqual(graph["a"]["c"]["weight"], 0.5)
        self.assertEqual(graph["b"]["c"]["weight"], 0.5)

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

    def test_hamming_distance_ignores_missing_values_by_default(self):
        from triplet_lineage.reconstruction import hamming_distance

        self.assertEqual(hamming_distance([1, -1, 0], [1, 2, 1]), 1)
        self.assertEqual(
            hamming_distance([1, -1, 0], [1, 2, 1], ignore_missing=False),
            2,
        )

    def test_ambiguous_triplet_is_abstained(self):
        import pandas as pd
        from triplet_lineage.reconstruction import infer_triplets_from_mutations

        class FakeCassiopeiaTree:
            leaves = ["a", "b", "c"]
            character_matrix = pd.DataFrame(
                [[1, -1], [1, -1], [1, -1]],
                index=leaves,
            )

        self.assertEqual(infer_triplets_from_mutations(FakeCassiopeiaTree()), [])

    def test_partition_triplets_recovers_consistent_top_level_split(self):
        from triplet_lineage.reconstruction import partition_triplets

        left, right = partition_triplets(
            [
                ("a", "b", "c"),
                ("a", "b", "d"),
                ("c", "d", "a"),
                ("c", "d", "b"),
            ],
            leaves=["a", "b", "c", "d"],
        )

        self.assertEqual({frozenset(left), frozenset(right)}, {frozenset({"a", "b"}), frozenset({"c", "d"})})

    def test_recursive_triplet_maxcut_tree_builds_hierarchy(self):
        import networkx as nx
        from triplet_lineage.reconstruction import recursive_triplet_maxcut_tree

        tree = recursive_triplet_maxcut_tree(
            leaves=["a", "b", "c", "d"],
            triplets=[
                ("a", "b", "c"),
                ("a", "b", "d"),
                ("c", "d", "a"),
                ("c", "d", "b"),
            ],
        )

        self.assertTrue(nx.is_arborescence(tree))
        self.assertEqual({node for node in tree.nodes if tree.out_degree(node) == 0}, {"a", "b", "c", "d"})
        cherry_sets = {
            frozenset(tree.successors(node))
            for node in tree.nodes
            if tree.out_degree(node) == 2
        }
        self.assertIn(frozenset({"a", "b"}), cherry_sets)
        self.assertIn(frozenset({"c", "d"}), cherry_sets)


if __name__ == "__main__":
    unittest.main()
