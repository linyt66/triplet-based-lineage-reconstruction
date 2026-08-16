"""Compare sampled Triplet Max-Cut with Cassiopeia baselines at larger scale."""

from __future__ import annotations

import argparse
import csv
import itertools
import time
from collections import defaultdict
from pathlib import Path

import cassiopeia.critique as critique
import cassiopeia.data as data
import networkx as nx
import numpy as np

from cassiopeia.solver.NeighborJoiningSolver import NeighborJoiningSolver
from cassiopeia.solver.UPGMASolver import UPGMASolver
from cassiopeia.solver.VanillaGreedySolver import VanillaGreedySolver

from triplet_lineage.metrics import calculate_triplets_correct
from triplet_lineage.simulation import (
    complete_binary_missing_tree_sim,
    complete_binary_tree_sim,
)

Q_DIST = {1: 0.98, 2: 0.01, 3: 0.01}
MISSING_STATE = -1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a 512-cell Cassiopeia benchmark for TMC and baselines."
    )
    parser.add_argument("--depth", type=int, default=10, help="Depth 10 gives 512 cells.")
    parser.add_argument("--sites", type=int, nargs="+", default=[120])
    parser.add_argument("--missing-rates", type=float, nargs="+", default=[0.1])
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--mutation-rate", type=float, default=0.8)
    parser.add_argument("--triplet-samples", type=int, default=300000)
    parser.add_argument("--accuracy-samples", type=int, default=20000)
    parser.add_argument("--tmc-restarts", type=int, default=8)
    parser.add_argument("--seed", type=int, default=5122026)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/cassiopeia_512cell_baseline_comparison.csv"),
    )
    return parser.parse_args()


def simulate_tree(k_sites: int, missing_rate: float, depth: int, mutation_rate: float):
    if missing_rate == 0:
        return complete_binary_tree_sim(
            k_cand=k_sites,
            q_dist=Q_DIST,
            lamb=mutation_rate,
            depth=depth,
        )
    return complete_binary_missing_tree_sim(
        k_cand=k_sites,
        q_dist=Q_DIST,
        lamb=mutation_rate,
        depth=depth,
        missing_rate=missing_rate,
    )


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


def build_signed_adjacency(triplets, subset: set[int]):
    adjacency = defaultdict(dict)
    for ingroup_a, ingroup_b, outgroup in triplets:
        if ingroup_a not in subset or ingroup_b not in subset or outgroup not in subset:
            continue
        add_edge(adjacency, ingroup_a, ingroup_b, -2.0)
        add_edge(adjacency, ingroup_a, outgroup, 1.0)
        add_edge(adjacency, ingroup_b, outgroup, 1.0)
    return adjacency


def add_edge(adjacency, left: int, right: int, weight: float) -> None:
    adjacency[left][right] = adjacency[left].get(right, 0.0) + weight
    adjacency[right][left] = adjacency[right].get(left, 0.0) + weight


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
                node_label = labels[node]
                for neighbor, weight in adjacency.get(node, {}).items():
                    if labels[neighbor] == node_label:
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


def balanced_split(nodes):
    ordered = sorted(nodes)
    midpoint = max(1, len(ordered) // 2)
    return set(ordered[:midpoint]), set(ordered[midpoint:])


def sampled_tmc_tree(leaves: list[str], triplets, rng: np.random.Generator, restarts: int):
    graph = nx.DiGraph()
    internal_id = itertools.count(1)

    def build(subset: set[int], node_name: str) -> str:
        graph.add_node(node_name)
        if len(subset) == 1:
            leaf = leaves[next(iter(subset))]
            graph.add_edge(node_name, leaf)
            return node_name
        if len(subset) == 2:
            for leaf_index in sorted(subset):
                graph.add_edge(node_name, leaves[leaf_index])
            return node_name

        adjacency = build_signed_adjacency(triplets, subset)
        left, right = local_search_maxcut(adjacency, sorted(subset), rng, restarts)
        if not left or not right:
            left, right = balanced_split(subset)

        for child_subset in (left, right):
            if len(child_subset) == 1:
                graph.add_edge(node_name, leaves[next(iter(child_subset))])
            else:
                child_name = f"tmc_internal_{next(internal_id)}"
                graph.add_edge(node_name, child_name)
                build(child_subset, child_name)
        return node_name

    build(set(range(len(leaves))), "root")
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


def normalized_rf_similarity(ground_tree, recon_tree) -> float:
    rf_distance, rf_max = critique.compare.robinson_foulds(ground_tree, recon_tree)
    if rf_max == 0:
        return 1.0
    return float(1 - rf_distance / rf_max)


def evaluate_networkx_tree(ground_tree, recon_graph, accuracy_samples: int) -> tuple[float, float]:
    triplet_accuracy = calculate_triplets_correct(ground_tree, recon_graph, n_samples=accuracy_samples)
    cass_tree = data.CassiopeiaTree(tree=recon_graph)
    return triplet_accuracy, normalized_rf_similarity(ground_tree, cass_tree)


def evaluate_cassiopeia_tree(ground_tree, recon_tree, accuracy_samples: int) -> tuple[float, float]:
    triplet_accuracy = calculate_triplets_correct(
        ground_tree,
        recon_tree.get_tree_topology(),
        n_samples=accuracy_samples,
    )
    return triplet_accuracy, normalized_rf_similarity(ground_tree, recon_tree)


def run_condition(args, k_sites: int, missing_rate: float, repeat: int) -> list[dict]:
    seed = args.seed + 100000 * repeat + 1000 * k_sites + int(100 * missing_rate)
    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    rows = []
    start = time.perf_counter()
    tree = simulate_tree(k_sites, missing_rate, args.depth, args.mutation_rate)
    simulation_seconds = time.perf_counter() - start
    leaves = list(tree.leaves)

    start = time.perf_counter()
    sampled_triplets = infer_sampled_triplets(tree, args.triplet_samples, rng)
    triplet_seconds = time.perf_counter() - start

    start = time.perf_counter()
    tmc_graph = sampled_tmc_tree(leaves, sampled_triplets, rng, args.tmc_restarts)
    tmc_solve_seconds = time.perf_counter() - start

    start = time.perf_counter()
    tmc_triplet_accuracy, tmc_rf_similarity = evaluate_networkx_tree(
        tree,
        tmc_graph,
        args.accuracy_samples,
    )
    tmc_eval_seconds = time.perf_counter() - start

    possible_triplets = len(leaves) * (len(leaves) - 1) * (len(leaves) - 2) // 6
    base = {
        "seed": seed,
        "depth": args.depth,
        "n_cells": len(leaves),
        "n_sites": k_sites,
        "missing_rate": missing_rate,
        "mutation_rate": args.mutation_rate,
        "repeat": repeat,
        "sampled_triplets_requested": args.triplet_samples,
        "sampled_triplets_inferred": len(sampled_triplets),
        "possible_triplets": possible_triplets,
        "sampled_triplet_coverage": len(sampled_triplets) / args.triplet_samples,
        "simulation_seconds": simulation_seconds,
    }
    rows.append(
        base
        | {
            "algorithm": "SampledRecursiveTMC",
            "triplet_accuracy": tmc_triplet_accuracy,
            "rf_similarity": tmc_rf_similarity,
            "triplet_inference_seconds": triplet_seconds,
            "solve_seconds": tmc_solve_seconds,
            "evaluation_seconds": tmc_eval_seconds,
            "status": "ok",
            "error": "",
        }
    )
    print(
        f"done SampledRecursiveTMC cells={len(leaves)} sites={k_sites} missing={missing_rate:.2f} "
        f"triplet_acc={tmc_triplet_accuracy:.3f} rf={tmc_rf_similarity:.3f} "
        f"sampled={len(sampled_triplets)} solve_sec={tmc_solve_seconds:.1f}",
        flush=True,
    )

    for algorithm in ("VanillaGreedy", "UPGMA", "NeighborJoining"):
        start = time.perf_counter()
        try:
            recon = solve_baseline(tree, algorithm)
            solve_seconds = time.perf_counter() - start
            eval_start = time.perf_counter()
            triplet_accuracy, rf_similarity = evaluate_cassiopeia_tree(
                tree,
                recon,
                args.accuracy_samples,
            )
            evaluation_seconds = time.perf_counter() - eval_start
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
            f"done {algorithm} cells={len(leaves)} sites={k_sites} missing={missing_rate:.2f} "
            f"triplet_acc={triplet_accuracy:.3f} rf={rf_similarity:.3f} "
            f"solve_sec={solve_seconds:.1f} status={status}",
            flush=True,
        )

    return rows


def main() -> None:
    args = parse_args()
    all_rows = []
    for k_sites in args.sites:
        for missing_rate in args.missing_rates:
            for repeat in range(args.repeats):
                all_rows.extend(run_condition(args, k_sites, missing_rate, repeat))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"wrote {args.output}")

    print("summary")
    for row in all_rows:
        print(
            f"{row['algorithm']:20s} cells={row['n_cells']} sites={row['n_sites']} "
            f"missing={row['missing_rate']:.2f} triplet={row['triplet_accuracy']:.3f} "
            f"rf={row['rf_similarity']:.3f} solve_sec={row['solve_seconds']:.1f} "
            f"status={row['status']}"
        )


if __name__ == "__main__":
    main()
