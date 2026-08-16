"""Stress-test sampled TMC on balanced 512-cell Cassiopeia simulations.

The stress settings target regimes where TMC should be most useful: low site
count, high dropout, high homoplasy, and conflicting local triplet evidence.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import time
from collections import defaultdict
from pathlib import Path

import cassiopeia.data as data
import cassiopeia.simulator as simulator
import networkx as nx
import numpy as np
from cassiopeia.solver import dissimilarity_functions
from cassiopeia.solver.NeighborJoiningSolver import NeighborJoiningSolver
from cassiopeia.solver.SharedMutationJoiningSolver import SharedMutationJoiningSolver
from cassiopeia.solver.UPGMASolver import UPGMASolver
from cassiopeia.solver.VanillaGreedySolver import VanillaGreedySolver

from triplet_lineage.metrics import calculate_triplets_correct

MISSING_STATE = -1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run balanced-topology stress benchmarks.")
    parser.add_argument("--n-cells", type=int, default=512)
    parser.add_argument("--sites", type=int, nargs="+", default=[20, 30])
    parser.add_argument("--edit-states", type=int, nargs="+", default=[3, 5])
    parser.add_argument("--missing-rates", type=float, nargs="+", default=[0.2, 0.3])
    parser.add_argument("--dominant-state-prob", type=float, default=0.0)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--mutation-rate", type=float, default=0.8)
    parser.add_argument("--triplet-samples", type=int, default=50000)
    parser.add_argument("--accuracy-samples", type=int, default=5000)
    parser.add_argument("--tmc-restarts", type=int, default=4)
    parser.add_argument("--refine-size", type=int, default=32)
    parser.add_argument("--min-signal-gap", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=805122026)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/cassiopeia_512cell_balanced_stress_benchmark.csv"),
    )
    return parser.parse_args()


def state_priors(n_edit_states: int, dominant_state_prob: float) -> dict[int, float]:
    if n_edit_states < 1:
        raise ValueError("n_edit_states must be positive.")
    if not 0 <= dominant_state_prob < 1:
        raise ValueError("dominant_state_prob must be in [0, 1).")
    if n_edit_states == 1:
        return {1: 1.0}
    if dominant_state_prob == 0:
        return {state: 1 / n_edit_states for state in range(1, n_edit_states + 1)}
    tail = (1 - dominant_state_prob) / (n_edit_states - 1)
    priors = {state: tail for state in range(1, n_edit_states + 1)}
    priors[1] = dominant_state_prob
    return priors


def balanced_random_topology(n_cells: int, rng: np.random.Generator) -> data.CassiopeiaTree:
    graph = nx.DiGraph()
    internal_id = 0
    leaf_id = 0

    def build(size: int) -> str:
        nonlocal internal_id, leaf_id
        if size == 1:
            node = f"c{leaf_id}"
            leaf_id += 1
            graph.add_node(node)
            return node
        node = "root" if internal_id == 0 else f"n{internal_id}"
        internal_id += 1
        low = max(1, int(round(size * 0.4)))
        high = min(size - 1, int(round(size * 0.6)))
        left_size = int(rng.integers(low, high + 1)) if low <= high else size // 2
        right_size = size - left_size
        left = build(left_size)
        right = build(right_size)
        graph.add_edge(node, left)
        graph.add_edge(node, right)
        return node

    build(n_cells)
    depths = nx.shortest_path_length(graph, source="root")
    leaves = {node for node in graph.nodes if graph.out_degree(node) == 0}
    internal_depths = [depth for node, depth in depths.items() if node not in leaves]
    max_internal_depth = max(internal_depths) if internal_depths else 1
    times = {}
    for node, depth in depths.items():
        if node == "root":
            times[node] = 0.0
        elif node in leaves:
            times[node] = 1.0
        else:
            times[node] = 0.95 * depth / max_internal_depth
    tree = data.CassiopeiaTree(tree=graph)
    tree.set_times(times)
    return tree


def overlay_cas9_data(
    tree,
    n_sites: int,
    n_edit_states: int,
    dominant_state_prob: float,
    mutation_rate: float,
    missing_rate: float,
):
    cas9 = simulator.Cas9LineageTracingDataSimulator(
        number_of_cassettes=n_sites,
        size_of_cassette=1,
        mutation_rate=mutation_rate,
        state_priors=state_priors(n_edit_states, dominant_state_prob),
        heritable_silencing_rate=0,
        stochastic_silencing_rate=missing_rate,
        heritable_missing_data_state=MISSING_STATE,
        stochastic_missing_data_state=MISSING_STATE,
    )
    cas9.overlay_data(tree)
    return tree


def observed_hamming_rate(values: np.ndarray, left: int, right: int) -> float:
    x = values[left]
    y = values[right]
    observed = (x != MISSING_STATE) & (y != MISSING_STATE)
    n_observed = int(np.sum(observed))
    if n_observed == 0:
        return np.inf
    return float(np.sum(x[observed] != y[observed]) / n_observed)


def random_distinct_triples(rng: np.random.Generator, n_items: int, n_samples: int) -> np.ndarray:
    triples = np.empty((0, 3), dtype=np.int32)
    while len(triples) < n_samples:
        batch_size = int((n_samples - len(triples)) * 1.3) + 1024
        batch = rng.integers(0, n_items, size=(batch_size, 3), dtype=np.int32)
        mask = (
            (batch[:, 0] != batch[:, 1])
            & (batch[:, 0] != batch[:, 2])
            & (batch[:, 1] != batch[:, 2])
        )
        triples = np.vstack([triples, batch[mask]])
    return triples[:n_samples]


def infer_confident_triplets(
    tree,
    n_samples: int,
    min_signal_gap: float,
    rng: np.random.Generator,
) -> list[tuple[int, int, int, float]]:
    leaves = list(tree.leaves)
    values = tree.character_matrix.loc[leaves].to_numpy()
    triples = random_distinct_triples(rng, len(leaves), n_samples)
    inferred = []
    for a, b, c in triples:
        distances = [
            (observed_hamming_rate(values, a, b), a, b, c),
            (observed_hamming_rate(values, a, c), a, c, b),
            (observed_hamming_rate(values, b, c), b, c, a),
        ]
        distances.sort(key=lambda item: item[0])
        best_distance, ingroup_a, ingroup_b, outgroup = distances[0]
        second_best = distances[1][0]
        gap = second_best - best_distance
        if np.isfinite(best_distance) and gap >= min_signal_gap:
            inferred.append((int(ingroup_a), int(ingroup_b), int(outgroup), float(gap)))
    return inferred


def add_edge(adjacency, left: int, right: int, weight: float) -> None:
    adjacency[left][right] = adjacency[left].get(right, 0.0) + weight
    adjacency[right][left] = adjacency[right].get(left, 0.0) + weight


def build_signed_adjacency(triplets, subset: set[int]):
    adjacency = defaultdict(dict)
    for ingroup_a, ingroup_b, outgroup, weight in triplets:
        if ingroup_a not in subset or ingroup_b not in subset or outgroup not in subset:
            continue
        add_edge(adjacency, ingroup_a, ingroup_b, -2.0 * weight)
        add_edge(adjacency, ingroup_a, outgroup, weight)
        add_edge(adjacency, ingroup_b, outgroup, weight)
    return adjacency


def cut_weight(adjacency, labels: dict[int, int]) -> float:
    score = 0.0
    seen = set()
    for left, neighbors in adjacency.items():
        for right, weight in neighbors.items():
            if (right, left) in seen:
                continue
            seen.add((left, right))
            if labels[left] != labels[right]:
                score += weight
    return score


def balanced_split(nodes):
    ordered = sorted(nodes)
    midpoint = max(1, len(ordered) // 2)
    return set(ordered[:midpoint]), set(ordered[midpoint:])


def local_search_maxcut(adjacency, nodes: list[int], rng: np.random.Generator, restarts: int):
    if len(nodes) < 2:
        return set(nodes), set()
    if not adjacency:
        return balanced_split(nodes)
    best_labels = None
    best_score = -np.inf
    ordered = sorted(nodes)
    for restart in range(restarts):
        shuffled = ordered.copy()
        rng.shuffle(shuffled)
        labels = {node: 1 if index < len(shuffled) // 2 else -1 for index, node in enumerate(shuffled)}
        if restart == 0:
            labels = {node: 1 if index < len(ordered) // 2 else -1 for index, node in enumerate(ordered)}
        improved = True
        while improved:
            improved = False
            best_node = None
            best_gain = 1e-12
            for node in ordered:
                gain = 0.0
                label = labels[node]
                for neighbor, weight in adjacency.get(node, {}).items():
                    gain += weight if labels[neighbor] == label else -weight
                if gain > best_gain:
                    best_gain = gain
                    best_node = node
            if best_node is not None:
                labels[best_node] *= -1
                improved = True
        score = cut_weight(adjacency, labels)
        if score > best_score:
            best_score = score
            best_labels = labels.copy()
    left = {node for node, label in best_labels.items() if label == 1}
    right = set(nodes) - left
    if not left or not right:
        return balanced_split(nodes)
    return left, right


def attach_smj_refinement(graph, parent, ground_tree, leaves, subset, prefix):
    if len(subset) == 1:
        graph.add_edge(parent, leaves[next(iter(subset))])
        return
    if len(subset) == 2:
        for leaf_index in sorted(subset):
            graph.add_edge(parent, leaves[leaf_index])
        return
    selected_leaves = [leaves[index] for index in sorted(subset)]
    sub_tree = data.CassiopeiaTree(
        character_matrix=ground_tree.character_matrix.loc[selected_leaves].astype(int),
        missing_state_indicator=MISSING_STATE,
    )
    solver = SharedMutationJoiningSolver(
        similarity_function=dissimilarity_functions.hamming_similarity_without_missing
    )
    solver.solve(sub_tree)
    topology = sub_tree.get_tree_topology()
    root = [node for node in topology.nodes if topology.in_degree(node) == 0][0]
    leaf_set = set(selected_leaves)

    def mapped(node):
        return node if node in leaf_set else f"{prefix}_smj_{node}"

    def copy_children(source_parent, target_parent):
        for child in topology.successors(source_parent):
            target_child = mapped(child)
            graph.add_edge(target_parent, target_child)
            if child not in leaf_set:
                copy_children(child, target_child)

    copy_children(root, parent)


def refined_tmc_tree(ground_tree, leaves, triplets, rng, restarts, refine_size):
    graph = nx.DiGraph()
    graph.add_node("root")
    internal_id = itertools.count(1)

    def build(parent, subset):
        if len(subset) <= refine_size:
            attach_smj_refinement(graph, parent, ground_tree, leaves, subset, parent)
            return
        adjacency = build_signed_adjacency(triplets, subset)
        left, right = local_search_maxcut(adjacency, sorted(subset), rng, restarts)
        for child_subset in (left, right):
            if len(child_subset) <= refine_size:
                attach_smj_refinement(
                    graph,
                    parent,
                    ground_tree,
                    leaves,
                    child_subset,
                    prefix=f"{parent}_{next(internal_id)}",
                )
            else:
                child_name = f"tmc_internal_{next(internal_id)}"
                graph.add_edge(parent, child_name)
                build(child_name, child_subset)

    build("root", set(range(len(leaves))))
    return graph


def fresh_tree_from_character_matrix(tree):
    return data.CassiopeiaTree(
        character_matrix=tree.character_matrix.astype(int),
        missing_state_indicator=MISSING_STATE,
    )


def solve_baseline(tree, algorithm: str):
    recon = fresh_tree_from_character_matrix(tree)
    if algorithm == "VanillaGreedy":
        solver = VanillaGreedySolver()
    elif algorithm == "UPGMA":
        solver = UPGMASolver(fast=False)
    elif algorithm == "NeighborJoining":
        solver = NeighborJoiningSolver(add_root=True, fast=False)
    else:
        raise ValueError(algorithm)
    solver.solve(recon)
    return recon


def rooted_clades(tree: nx.DiGraph, leaves: set[str]) -> set[frozenset[str]]:
    clades = set()
    for node in tree.nodes:
        descendants = nx.descendants(tree, node)
        leaf_descendants = frozenset(descendants & leaves)
        if 1 < len(leaf_descendants) < len(leaves):
            clades.add(leaf_descendants)
    return clades


def rooted_rf_similarity(ground_graph: nx.DiGraph, recon_graph: nx.DiGraph) -> float:
    ground_leaves = {node for node in ground_graph.nodes if ground_graph.out_degree(node) == 0}
    recon_leaves = {node for node in recon_graph.nodes if recon_graph.out_degree(node) == 0}
    leaves = ground_leaves & recon_leaves
    if not leaves:
        return 0.0
    ground_clades = rooted_clades(ground_graph, leaves)
    recon_clades = rooted_clades(recon_graph, leaves)
    denominator = len(ground_clades) + len(recon_clades)
    if denominator == 0:
        return 1.0
    rf_distance = len(ground_clades - recon_clades) + len(recon_clades - ground_clades)
    return float(1 - rf_distance / denominator)


def evaluate_networkx_tree(ground_tree, recon_graph, accuracy_samples: int):
    triplet_accuracy = calculate_triplets_correct(ground_tree, recon_graph, n_samples=accuracy_samples)
    return triplet_accuracy, rooted_rf_similarity(ground_tree.get_tree_topology(), recon_graph)


def evaluate_cassiopeia_tree(ground_tree, recon_tree, accuracy_samples: int):
    recon_graph = recon_tree.get_tree_topology()
    triplet_accuracy = calculate_triplets_correct(ground_tree, recon_graph, n_samples=accuracy_samples)
    return triplet_accuracy, rooted_rf_similarity(ground_tree.get_tree_topology(), recon_graph)


def run_condition(args, condition_index, n_sites, n_edit_states, missing_rate, repeat):
    seed = args.seed + 100000 * condition_index + repeat
    rng = np.random.default_rng(seed)
    np.random.seed(seed)
    rows = []

    start = time.perf_counter()
    tree = balanced_random_topology(args.n_cells, rng)
    overlay_cas9_data(
        tree,
        n_sites,
        n_edit_states,
        args.dominant_state_prob,
        args.mutation_rate,
        missing_rate,
    )
    simulation_seconds = time.perf_counter() - start
    leaves = list(tree.leaves)

    start = time.perf_counter()
    triplets = infer_confident_triplets(tree, args.triplet_samples, args.min_signal_gap, rng)
    triplet_seconds = time.perf_counter() - start
    possible_triplets = len(leaves) * (len(leaves) - 1) * (len(leaves) - 2) // 6
    base = {
        "seed": seed,
        "topology_model": "balanced_random_40_60_splits",
        "n_cells": len(leaves),
        "n_sites": n_sites,
        "edit_states": n_edit_states,
        "dominant_state_prob": args.dominant_state_prob,
        "missing_rate": missing_rate,
        "mutation_rate": args.mutation_rate,
        "repeat": repeat,
        "triplet_samples_requested": args.triplet_samples,
        "sampled_triplets_inferred": len(triplets),
        "possible_triplets": possible_triplets,
        "sampled_triplet_coverage": len(triplets) / args.triplet_samples,
        "min_signal_gap": args.min_signal_gap,
        "simulation_seconds": simulation_seconds,
    }

    start = time.perf_counter()
    tmc_graph = refined_tmc_tree(tree, leaves, triplets, rng, args.tmc_restarts, args.refine_size)
    solve_seconds = time.perf_counter() - start
    triplet_accuracy, rf_similarity = evaluate_networkx_tree(tree, tmc_graph, args.accuracy_samples)
    rows.append(base | {
        "algorithm": f"ConfidentWeightedTMC_SMJ{args.refine_size}",
        "triplet_accuracy": triplet_accuracy,
        "rf_similarity": rf_similarity,
        "triplet_inference_seconds": triplet_seconds,
        "solve_seconds": solve_seconds,
        "status": "ok",
        "error": "",
    })
    print(
        f"condition={condition_index} repeat={repeat} TMC sites={n_sites} states={n_edit_states} "
        f"missing={missing_rate:.2f} coverage={base['sampled_triplet_coverage']:.3f} "
        f"triplet={triplet_accuracy:.3f} rf={rf_similarity:.3f} solve={solve_seconds:.1f}s",
        flush=True,
    )

    for algorithm in ("VanillaGreedy", "UPGMA", "NeighborJoining"):
        start = time.perf_counter()
        try:
            recon = solve_baseline(tree, algorithm)
            baseline_solve = time.perf_counter() - start
            triplet_accuracy, rf_similarity = evaluate_cassiopeia_tree(
                tree,
                recon,
                args.accuracy_samples,
            )
            status = "ok"
            error = ""
        except Exception as exc:
            baseline_solve = time.perf_counter() - start
            triplet_accuracy = np.nan
            rf_similarity = np.nan
            status = "error"
            error = f"{type(exc).__name__}: {exc}"
        rows.append(base | {
            "algorithm": algorithm,
            "triplet_accuracy": triplet_accuracy,
            "rf_similarity": rf_similarity,
            "triplet_inference_seconds": 0.0,
            "solve_seconds": baseline_solve,
            "status": status,
            "error": error,
        })
        print(
            f"condition={condition_index} repeat={repeat} {algorithm} sites={n_sites} states={n_edit_states} "
            f"missing={missing_rate:.2f} triplet={triplet_accuracy:.3f} rf={rf_similarity:.3f} "
            f"solve={baseline_solve:.1f}s status={status}",
            flush=True,
        )
    return rows


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    all_rows = []
    condition_index = 0
    for n_sites in args.sites:
        for n_edit_states in args.edit_states:
            for missing_rate in args.missing_rates:
                for repeat in range(args.repeats):
                    all_rows.extend(
                        run_condition(args, condition_index, n_sites, n_edit_states, missing_rate, repeat)
                    )
                    write_rows(args.output, all_rows)
                condition_index += 1
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
