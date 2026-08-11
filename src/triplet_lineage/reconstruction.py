"""Triplet-based and baseline lineage reconstruction routines."""

from typing import Iterable, List, Sequence, Set, Tuple

import cassiopeia.data as data
import networkx as nx
import numpy as np
from cassiopeia.solver import dissimilarity_functions
from cassiopeia.solver.PercolationSolver import PercolationSolver
from cassiopeia.solver.SharedMutationJoiningSolver import SharedMutationJoiningSolver
from cassiopeia.solver.VanillaGreedySolver import VanillaGreedySolver

import cvxgraphalgs as cvxgr

Triplet = Tuple[str, str, str]


def percolation_solve(tree: data.CassiopeiaTree) -> data.CassiopeiaTree:
    """Reconstruct a tree using Cassiopeia's percolation solver."""
    cm = tree.character_matrix.astype(int)
    recon_tree = data.CassiopeiaTree(
        character_matrix=cm,
        missing_state_indicator=-1,
    )

    solver = PercolationSolver(joining_solver=VanillaGreedySolver())
    solver.solve(recon_tree)
    return recon_tree


def shared_mutation_solve(tree: data.CassiopeiaTree) -> data.CassiopeiaTree:
    """Reconstruct a tree using Shared Mutation Joining."""
    cm = tree.character_matrix.astype(int)
    recon_tree = data.CassiopeiaTree(
        character_matrix=cm,
        missing_state_indicator=-1,
    )

    solver = SharedMutationJoiningSolver(
        similarity_function=dissimilarity_functions.hamming_similarity_without_missing
    )
    solver.solve(recon_tree)
    return recon_tree


def hamming_distance(x: Sequence[int], y: Sequence[int]) -> int:
    """Compute the Hamming distance between two mutation vectors."""
    return int(np.sum(np.asarray(x) != np.asarray(y)))


def infer_triplets_from_mutations(tree: data.CassiopeiaTree) -> List[Triplet]:
    """Infer rooted triplets from the closest pair in each three-leaf set."""
    cm = tree.character_matrix
    triplets = []

    for a, b, c in _leaf_triplets(tree.leaves):
        distances = {
            (a, b): hamming_distance(cm.loc[a], cm.loc[b]),
            (a, c): hamming_distance(cm.loc[a], cm.loc[c]),
            (b, c): hamming_distance(cm.loc[b], cm.loc[c]),
        }
        (ingroup_a, ingroup_b), _ = min(distances.items(), key=lambda item: item[1])
        outgroup = ({a, b, c} - {ingroup_a, ingroup_b}).pop()
        triplets.append((ingroup_a, ingroup_b, outgroup))

    return triplets


def construct_triplet_graph(triplets: Iterable[Triplet]) -> nx.Graph:
    """Encode triplet constraints as a weighted graph for Max-Cut."""
    graph = nx.Graph()

    for ingroup_a, ingroup_b, outgroup in triplets:
        _add_weighted_edge(graph, ingroup_a, ingroup_b, -2)
        _add_weighted_edge(graph, ingroup_a, outgroup, 1)
        _add_weighted_edge(graph, ingroup_b, outgroup, 1)

    return graph


def gw_partition(triplets: Iterable[Triplet]) -> Tuple[Set[str], Set[str]]:
    """Partition leaves with the Goemans-Williamson SDP relaxation."""
    graph = construct_triplet_graph(triplets)
    cut = cvxgr.algorithms.goemans_williamson_weighted(graph)
    return set(cut.left), set(cut.right)


def goemans_williamson_solve_desired_triplets(
    triplets: Iterable[Triplet],
) -> Tuple[Set[str], Set[str]]:
    """Compatibility wrapper for the original notebook function name."""
    return gw_partition(triplets)


def build_tree_from_triplets(
    tree: data.CassiopeiaTree,
    triplets: Iterable[Triplet],
) -> nx.DiGraph:
    """Build the TMC reconstruction from a triplet partition."""
    left, right = gw_partition(triplets)
    cm = tree.character_matrix

    subtrees = []
    for subset in (left, right):
        sub_cm = cm.loc[list(subset)]
        sub_tree = data.CassiopeiaTree(
            character_matrix=sub_cm,
            missing_state_indicator=-1,
        )
        subtrees.append(shared_mutation_solve(sub_tree).get_tree_topology())

    recon_tree = nx.DiGraph()
    root = "root"
    recon_tree.add_node(root)

    for subtree in subtrees:
        recon_tree.add_nodes_from(subtree.nodes)
        recon_tree.add_edges_from(subtree.edges)
        recon_tree.add_edge(root, list(subtree.nodes)[0])

    return recon_tree


def build_tree_from_triplet_partition(
    tree: data.CassiopeiaTree,
    triplets: Iterable[Triplet],
) -> nx.DiGraph:
    """Compatibility wrapper for the original TMC reconstruction name."""
    return build_tree_from_triplets(tree, triplets)


def _add_weighted_edge(graph: nx.Graph, u: str, v: str, delta: int) -> None:
    current = graph.get_edge_data(u, v, {"weight": 0})["weight"]
    graph.add_edge(u, v, weight=current + delta)


def _leaf_triplets(leaves):
    import itertools

    return itertools.combinations(leaves, 3)


find_recon_triplets = infer_triplets_from_mutations
