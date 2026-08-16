"""Run 512-cell, 30-site repeated benchmarks on random topologies.

The output is a long-format CSV with one row per repeat and algorithm, suitable
for violin plots of triplet accuracy or RF similarity.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import time
from collections import defaultdict
from pathlib import Path

import cassiopeia.critique as critique
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
    parser = argparse.ArgumentParser(
        description="Benchmark 512-cell, 30-site Cassiopeia data on random topologies."
    )
    parser.add_argument("--n-cells", type=int, default=512)
    parser.add_argument("--sites", type=int, default=30)
    parser.add_argument("--missing-rate", type=float, default=0.1)
    parser.add_argument("--edit-states", type=int, default=30)
    parser.add_argument("--repeats", type=int, default=50)
    parser.add_argument("--mutation-rate", type=float, default=0.8)
    parser.add_argument("--triplet-samples", type=int, default=100000)
    parser.add_argument("--accuracy-samples", type=int, default=10000)
    parser.add_argument("--tmc-restarts", type=int, default=6)
    parser.add_argument("--refine-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=305122026)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/cassiopeia_512cell_30site_random_topology_50repeats.csv"),
    )
    return parser.parse_args()


def random_binary_topology(n_cells: int, rng: np.random.Generator) -> data.CassiopeiaTree:
    graph = nx.DiGraph()
    root = "root"
    graph.add_node(root)
    active = [root]
    internal_id = itertools.count(1)

    while len(active) < n_cells:
        parent_index = int(rng.integers(0, len(active)))
        parent = active.pop(parent_index)
        left = f"n{next(internal_id)}"
        right = f"n{next(internal_id)}"
        graph.add_edge(parent, left)
        graph.add_edge(parent, right)
        active.extend([left, right])

    mapping = {node: f"c{index}" for index, node in enumerate(sorted(active, key=str))}
    graph = nx.relabel_nodes(graph, mapping, copy=True)
    leaves = set(mapping.values())

    depths = nx.shortest_path_length(graph, source=root)
    internal_depths = [depth for node, depth in depths.items() if node not in leaves]
    max_internal_depth = max(internal_depths) if internal_depths else 1
    times = {}
    for node, depth in depths.items():
        if node == root:
            times[node] = 0.0
        elif node in leaves:
            times[node] = 1.0
        else:
            times[node] = 0.95 * depth / max_internal_depth

    tree = data.CassiopeiaTree(tree=graph)
    tree.set_times(times)
    return tree


def uniform_state_priors(n_edit_states: int) -> dict[int, float]:
    if n_edit_states < 1:
        raise ValueError("n_edit_states must be positive.")
    return {state: 1 / n_edit_states for state in range(1, n_edit_states + 1)}


def overlay_cas9_data(
    tree,
    n_sites: int,
    mutation_rate: float,
    missing_rate: float,
    n_edit_states: int,
):
    cas9 = simulator.Cas9LineageTracingDataSimulator(
        number_of_cassettes=n_sites,
        size_of_cassette=1,
        mutation_rate=mutation_rate,
        state_priors=uniform_state_priors(n_edit_states),
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


def infer_sampled_triplets(tree, n_samples: int, rng: np.random.Generator) -> list[tuple[int, int, int]]:
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
        if np.isfinite(best_distance) and second_best > best_distance:
            inferred.append((int(ingroup_a), int(ingroup_b), int(outgroup)))

    return inferred


def add_edge(adjacency, left: int, right: int, weight: float) -> None:
    adjacency[left][right] = adjacency[left].get(right, 0.0) + weight
    adjacency[right][left] = adjacency[right].get(left, 0.0) + weight


def build_signed_adjacency(triplets, subset: set[int]):
    adjacency = defaultdict(dict)
    for ingroup_a, ingroup_b, outgroup in triplets:
        if ingroup_a not in subset or ingroup_b not in subset or outgroup not in subset:
            continue
        add_edge(adjacency, ingroup_a, ingroup_b, -2.0)
        add_edge(adjacency, ingroup_a, outgroup, 1.0)
        add_edge(adjacency, ingroup_b, outgroup, 1.0)
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
        labels = {}
        midpoint = len(shuffled) // 2
        for index, node in enumerate(shuffled):
            labels[node] = 1 if index < midpoint else -1
        if restart == 0:
            for index, node in enumerate(ordered):
                labels[node] = 1 if index < len(ordered) // 2 else -1

        improved = True
        while improved:
            improved = False
            best_node = None
            best_gain = 1e-12
            for node in ordered:
                gain = 0.0
                label = labels[node]
                for neighbor, weight in adjacency.get(node, {}).items():
                    if labels[neighbor] == label:
                        gain += weight
                    else:
                        gain -= weight
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
    roots = [node for node in topology.nodes if topology.in_degree(node) == 0]
    sub_root = roots[0]
    leaf_set = set(selected_leaves)

    def mapped(node):
        if node in leaf_set:
            return node
        return f"{prefix}_smj_{node}"

    def copy_children(source_parent, target_parent):
        for child in topology.successors(source_parent):
            target_child = mapped(child)
            graph.add_edge(target_parent, target_child)
            if child not in leaf_set:
                copy_children(child, target_child)

    copy_children(sub_root, parent)


def refined_sampled_tmc_tree(ground_tree, leaves, triplets, rng, restarts, refine_size):
    graph = nx.DiGraph()
    graph.add_node("root")
    internal_id = itertools.count(1)

    def build(parent, subset):
        if len(subset) <= refine_size:
            attach_smj_refinement(graph, parent, ground_tree, leaves, subset, parent)
            return

        adjacency = build_signed_adjacency(triplets, subset)
        left, right = local_search_maxcut(adjacency, sorted(subset), rng, restarts)
        if not left or not right:
            left, right = balanced_split(subset)

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
    """Return non-trivial rooted leaf clades for a directed tree."""
    clades = set()
    for node in tree.nodes:
        descendants = nx.descendants(tree, node)
        leaf_descendants = frozenset(descendants & leaves)
        if 1 < len(leaf_descendants) < len(leaves):
            clades.add(leaf_descendants)
    return clades


def rooted_rf_similarity(ground_graph: nx.DiGraph, recon_graph: nx.DiGraph) -> float:
    """Compute rooted RF similarity from explicit clade sets."""
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
    triplet_accuracy = calculate_triplets_correct(
        ground_tree,
        recon_graph,
        n_samples=accuracy_samples,
    )
    return triplet_accuracy, rooted_rf_similarity(ground_tree.get_tree_topology(), recon_graph)


def run_repeat(args, repeat: int) -> list[dict]:
    seed = args.seed + repeat
    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    start = time.perf_counter()
    tree = random_binary_topology(args.n_cells, rng)
    overlay_cas9_data(
        tree,
        args.sites,
        args.mutation_rate,
        args.missing_rate,
        args.edit_states,
    )
    simulation_seconds = time.perf_counter() - start

    leaves = list(tree.leaves)
    possible_triplets = len(leaves) * (len(leaves) - 1) * (len(leaves) - 2) // 6
    base = {
        "seed": seed,
        "topology_model": "random_binary_leaf_splitting",
        "n_cells": len(leaves),
        "n_sites": args.sites,
        "missing_rate": args.missing_rate,
        "mutation_rate": args.mutation_rate,
        "edit_states": args.edit_states,
        "repeat": repeat,
        "triplet_samples_requested": args.triplet_samples,
        "possible_triplets": possible_triplets,
        "simulation_seconds": simulation_seconds,
    }

    rows = []
    start = time.perf_counter()
    sampled_triplets = infer_sampled_triplets(tree, args.triplet_samples, rng)
    triplet_inference_seconds = time.perf_counter() - start
    base["sampled_triplets_inferred"] = len(sampled_triplets)
    base["sampled_triplet_coverage"] = len(sampled_triplets) / args.triplet_samples

    start = time.perf_counter()
    tmc_graph = refined_sampled_tmc_tree(
        tree,
        leaves,
        sampled_triplets,
        rng,
        args.tmc_restarts,
        args.refine_size,
    )
    solve_seconds = time.perf_counter() - start
    start = time.perf_counter()
    triplet_accuracy, rf_similarity = evaluate_networkx_tree(
        tree,
        tmc_graph,
        args.accuracy_samples,
    )
    evaluation_seconds = time.perf_counter() - start
    rows.append(
        base
        | {
            "algorithm": f"SampledRecursiveTMC_SMJ{args.refine_size}",
            "triplet_accuracy": triplet_accuracy,
            "rf_similarity": rf_similarity,
            "triplet_inference_seconds": triplet_inference_seconds,
            "solve_seconds": solve_seconds,
            "evaluation_seconds": evaluation_seconds,
            "status": "ok",
            "error": "",
        }
    )
    print(
        f"repeat={repeat:03d} SampledRecursiveTMC_SMJ{args.refine_size} "
        f"triplet={triplet_accuracy:.3f} rf={rf_similarity:.3f} "
        f"solve={solve_seconds:.1f}s",
        flush=True,
    )

    for algorithm in ("VanillaGreedy", "UPGMA", "NeighborJoining"):
        start = time.perf_counter()
        try:
            recon = solve_baseline(tree, algorithm)
            solve_seconds = time.perf_counter() - start
            start = time.perf_counter()
            triplet_accuracy, rf_similarity = evaluate_cassiopeia_tree(
                tree,
                recon,
                args.accuracy_samples,
            )
            evaluation_seconds = time.perf_counter() - start
            status = "ok"
            error = ""
        except Exception as exc:
            solve_seconds = time.perf_counter() - start
            evaluation_seconds = 0.0
            triplet_accuracy = np.nan
            rf_similarity = np.nan
            status = "error"
            error = f"{type(exc).__name__}: {exc}"
        rows.append(
            base
            | {
                "algorithm": algorithm,
                "triplet_accuracy": triplet_accuracy,
                "rf_similarity": rf_similarity,
                "triplet_inference_seconds": 0.0,
                "solve_seconds": solve_seconds,
                "evaluation_seconds": evaluation_seconds,
                "status": status,
                "error": error,
            }
        )
        print(
            f"repeat={repeat:03d} {algorithm} triplet={triplet_accuracy:.3f} "
            f"rf={rf_similarity:.3f} solve={solve_seconds:.1f}s status={status}",
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
    started = time.perf_counter()
    for repeat in range(args.repeats):
        all_rows.extend(run_repeat(args, repeat))
        write_rows(args.output, all_rows)
        print(f"checkpoint wrote {len(all_rows)} rows to {args.output}", flush=True)

    elapsed = time.perf_counter() - started
    print(f"completed {args.repeats} repeats in {elapsed:.1f}s")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()





