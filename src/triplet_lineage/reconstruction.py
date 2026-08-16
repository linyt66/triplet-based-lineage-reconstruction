"""Triplet-based and baseline lineage reconstruction routines."""

from __future__ import annotations

from itertools import count
import math
from typing import TYPE_CHECKING, Iterable, List, Sequence, Set, Tuple

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    import cassiopeia.data as data

Triplet = Tuple[str, str, str]
WeightedTriplet = Tuple[str, str, str, float]
TripletLike = Triplet | WeightedTriplet


def percolation_solve(tree: data.CassiopeiaTree) -> data.CassiopeiaTree:
    """Reconstruct a tree using Cassiopeia's percolation solver."""
    import cassiopeia.data as data
    from cassiopeia.solver.PercolationSolver import PercolationSolver
    from cassiopeia.solver.VanillaGreedySolver import VanillaGreedySolver

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
    import cassiopeia.data as data
    from cassiopeia.solver import dissimilarity_functions
    from cassiopeia.solver.SharedMutationJoiningSolver import SharedMutationJoiningSolver

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


def hamming_distance(
    x: Sequence[int],
    y: Sequence[int],
    missing_state: int | None = -1,
    ignore_missing: bool = True,
) -> int:
    """Compute Hamming distance between two mutation vectors.

    Missing values are ignored by default so that dropout does not create a
    false disagreement between otherwise identical barcodes.
    """
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)

    if missing_state is not None and ignore_missing:
        observed = (x_arr != missing_state) & (y_arr != missing_state)
        x_arr = x_arr[observed]
        y_arr = y_arr[observed]

    return int(np.sum(x_arr != y_arr))


def observed_hamming_rate(
    x: Sequence[int],
    y: Sequence[int],
    missing_state: int = -1,
) -> float:
    """Return normalized Hamming distance on jointly observed sites."""
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)
    observed = (x_arr != missing_state) & (y_arr != missing_state)
    n_observed = int(np.sum(observed))
    if n_observed == 0:
        return math.inf
    return float(np.sum(x_arr[observed] != y_arr[observed]) / n_observed)


def infer_triplets_from_mutations(
    tree: data.CassiopeiaTree,
    missing_state: int = -1,
    min_signal_gap: float = 0.0,
) -> List[Triplet]:
    """Infer rooted triplets from the closest observed barcode pair.

    Ties and weakly separated triplets are treated as abstentions. This mirrors
    the manuscript's distinction between resolvable triplets and triplets whose
    local signal gap is too small to support a confident rooted statement.
    """
    cm = tree.character_matrix
    triplets = []

    for a, b, c in _leaf_triplets(tree.leaves):
        distances = {
            (a, b): observed_hamming_rate(cm.loc[a], cm.loc[b], missing_state),
            (a, c): observed_hamming_rate(cm.loc[a], cm.loc[c], missing_state),
            (b, c): observed_hamming_rate(cm.loc[b], cm.loc[c], missing_state),
        }
        ranked = sorted(distances.items(), key=lambda item: item[1])
        (ingroup_a, ingroup_b), best_distance = ranked[0]
        second_best_distance = ranked[1][1]
        if not math.isfinite(best_distance):
            continue
        if second_best_distance - best_distance <= min_signal_gap:
            continue

        outgroup = ({a, b, c} - {ingroup_a, ingroup_b}).pop()
        triplets.append((ingroup_a, ingroup_b, outgroup))

    return triplets


def construct_triplet_graph(triplets: Iterable[TripletLike]) -> nx.Graph:
    """Encode triplet constraints as a weighted graph for Max-Cut.

    A triplet ``(a, b, c)`` means ``a`` and ``b`` should remain on the same side
    of the cut, while ``c`` should be separated. An optional fourth value is
    interpreted as a confidence weight.
    """
    graph = nx.Graph()

    for triplet in triplets:
        ingroup_a, ingroup_b, outgroup, weight = _unpack_triplet(triplet)
        _add_weighted_edge(graph, ingroup_a, ingroup_b, -2 * weight)
        _add_weighted_edge(graph, ingroup_a, outgroup, weight)
        _add_weighted_edge(graph, ingroup_b, outgroup, weight)

    return graph


def exact_maxcut_partition(graph: nx.Graph) -> Tuple[Set[str], Set[str]]:
    """Find the exact maximum cut for a small weighted graph."""
    nodes = sorted(graph.nodes, key=str)
    if len(nodes) < 2:
        return set(nodes), set()

    best_left: Set[str] | None = None
    best_score = -math.inf
    best_balance = -math.inf
    n_nodes = len(nodes)

    # Fix the first node on the left to remove the symmetric duplicate cut.
    for mask in range(1, 2 ** (n_nodes - 1)):
        left = {nodes[0]}
        right = set()
        for offset, node in enumerate(nodes[1:]):
            if mask & (1 << offset):
                right.add(node)
            else:
                left.add(node)

        score = _cut_weight(graph, left, right)
        balance = -abs(len(left) - len(right))
        if score > best_score or (score == best_score and balance > best_balance):
            best_score = score
            best_balance = balance
            best_left = left

    if best_left is None:
        return {nodes[0]}, set(nodes[1:])

    return best_left, set(nodes) - best_left


def partition_triplets(
    triplets: Iterable[TripletLike],
    leaves: Iterable[str] | None = None,
    exact_threshold: int = 12,
) -> Tuple[Set[str], Set[str]]:
    """Partition leaves by maximizing triplet-consistency cut weight."""
    triplet_list = list(triplets)
    graph = construct_triplet_graph(triplet_list)
    if leaves is not None:
        graph.add_nodes_from(leaves)

    nodes = sorted(graph.nodes, key=str)
    if len(nodes) < 2:
        return set(nodes), set()
    if graph.number_of_edges() == 0:
        return _balanced_split(nodes)
    if len(nodes) <= exact_threshold:
        return exact_maxcut_partition(graph)

    try:
        return gw_partition_from_graph(graph)
    except Exception:
        return _balanced_split(nodes)


def gw_partition(triplets: Iterable[TripletLike]) -> Tuple[Set[str], Set[str]]:
    """Partition leaves with the Goemans-Williamson SDP relaxation."""
    graph = construct_triplet_graph(triplets)
    return gw_partition_from_graph(graph)


def gw_partition_from_graph(graph: nx.Graph) -> Tuple[Set[str], Set[str]]:
    """Partition a weighted graph with the Goemans-Williamson relaxation."""
    import cvxgraphalgs as cvxgr

    cut = cvxgr.algorithms.goemans_williamson_weighted(graph)
    return set(cut.left), set(cut.right)


def goemans_williamson_solve_desired_triplets(
    triplets: Iterable[TripletLike],
) -> Tuple[Set[str], Set[str]]:
    """Compatibility wrapper for the original notebook function name."""
    return gw_partition(triplets)


def build_tree_from_triplets(
    tree: data.CassiopeiaTree,
    triplets: Iterable[TripletLike],
    exact_threshold: int = 12,
) -> nx.DiGraph:
    """Build a recursive Triplet Max-Cut reconstruction.

    The manuscript describes TMC as a top-down recursive Max-Cut procedure. This
    implementation follows that description directly, filtering the triplet set
    to the active leaf subset at each recursive split.
    """
    return recursive_triplet_maxcut_tree(
        leaves=tree.leaves,
        triplets=triplets,
        exact_threshold=exact_threshold,
    )


def recursive_triplet_maxcut_tree(
    leaves: Iterable[str],
    triplets: Iterable[TripletLike],
    exact_threshold: int = 12,
    root_name: str = "root",
) -> nx.DiGraph:
    """Construct a rooted tree by recursively applying Triplet Max-Cut."""
    leaf_list = list(leaves)
    triplet_list = list(triplets)
    leaf_set = set(leaf_list)
    root = _unique_name(root_name, leaf_set)
    internal_ids = count(1)
    graph = nx.DiGraph()

    def build(subset: Set[str], node_name: str | None = None) -> str:
        if len(subset) == 1:
            leaf = next(iter(subset))
            graph.add_node(leaf)
            return leaf

        if node_name is None:
            node_name = _unique_name(f"tmc_internal_{next(internal_ids)}", leaf_set)
        graph.add_node(node_name)

        if len(subset) == 2:
            for leaf in sorted(subset, key=str):
                graph.add_edge(node_name, leaf)
            return node_name

        local_triplets = _triplets_within_subset(triplet_list, subset)
        left, right = partition_triplets(
            local_triplets,
            leaves=subset,
            exact_threshold=exact_threshold,
        )
        if not left or not right:
            left, right = _balanced_split(sorted(subset, key=str))

        for child_subset in (left, right):
            child = build(set(child_subset))
            graph.add_edge(node_name, child)
        return node_name

    build(leaf_set, root)
    return graph


def single_level_maxcut_smj_tree(
    tree: data.CassiopeiaTree,
    triplets: Iterable[TripletLike],
    exact_threshold: int = 12,
) -> nx.DiGraph:
    """Build the single-level Max-Cut plus SMJ ablation tree.

    This helper preserves the earlier non-recursive baseline so validation
    notebooks can isolate the gain from recursive TMC.
    """
    import cassiopeia.data as data

    left, right = partition_triplets(
        triplets,
        leaves=tree.leaves,
        exact_threshold=exact_threshold,
    )
    cm = tree.character_matrix

    recon_tree = nx.DiGraph()
    root = _unique_name("root", set(tree.leaves))
    recon_tree.add_node(root)

    for subset in (left, right):
        if not subset:
            continue
        if len(subset) == 1:
            leaf = next(iter(subset))
            recon_tree.add_edge(root, leaf)
            continue

        sub_cm = cm.loc[sorted(subset, key=str)]
        sub_tree = data.CassiopeiaTree(
            character_matrix=sub_cm,
            missing_state_indicator=-1,
        )
        subtree = shared_mutation_solve(sub_tree).get_tree_topology()
        recon_tree.add_nodes_from(subtree.nodes)
        recon_tree.add_edges_from(subtree.edges)
        subtree_root = next(
            node for node in subtree.nodes if subtree.in_degree(node) == 0
        )
        recon_tree.add_edge(root, subtree_root)

    return recon_tree
def build_tree_from_triplet_partition(
    tree: data.CassiopeiaTree,
    triplets: Iterable[TripletLike],
) -> nx.DiGraph:
    """Compatibility wrapper for the original TMC reconstruction name."""
    return build_tree_from_triplets(tree, triplets)


def _add_weighted_edge(graph: nx.Graph, u: str, v: str, delta: float) -> None:
    current = graph.get_edge_data(u, v, {"weight": 0})["weight"]
    graph.add_edge(u, v, weight=current + delta)


def _unpack_triplet(triplet: TripletLike) -> WeightedTriplet:
    if len(triplet) == 3:
        ingroup_a, ingroup_b, outgroup = triplet
        return ingroup_a, ingroup_b, outgroup, 1.0
    if len(triplet) == 4:
        ingroup_a, ingroup_b, outgroup, weight = triplet
        return ingroup_a, ingroup_b, outgroup, float(weight)
    raise ValueError("Triplets must have three labels and an optional weight.")


def _cut_weight(graph: nx.Graph, left: Set[str], right: Set[str]) -> float:
    return float(
        sum(
            data.get("weight", 1.0)
            for u, v, data in graph.edges(data=True)
            if (u in left and v in right) or (u in right and v in left)
        )
    )


def _balanced_split(nodes: Iterable[str]) -> Tuple[Set[str], Set[str]]:
    ordered = sorted(nodes, key=str)
    midpoint = max(1, len(ordered) // 2)
    return set(ordered[:midpoint]), set(ordered[midpoint:])


def _triplets_within_subset(
    triplets: Iterable[TripletLike],
    subset: Set[str],
) -> List[TripletLike]:
    return [
        triplet
        for triplet in triplets
        if set(_unpack_triplet(triplet)[:3]).issubset(subset)
    ]


def _unique_name(base: str, unavailable: Set[str]) -> str:
    if base not in unavailable:
        return base

    suffix = 1
    candidate = f"{base}_{suffix}"
    while candidate in unavailable:
        suffix += 1
        candidate = f"{base}_{suffix}"
    return candidate


def _leaf_triplets(leaves):
    import itertools

    return itertools.combinations(leaves, 3)


find_recon_triplets = infer_triplets_from_mutations

