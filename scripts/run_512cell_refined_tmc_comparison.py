"""Run the 512-cell benchmark with sampled TMC and SMJ local refinement."""

from __future__ import annotations

import argparse
import csv
import itertools
import time
from pathlib import Path

import cassiopeia.data as data
import networkx as nx
import numpy as np
from cassiopeia.solver import dissimilarity_functions
from cassiopeia.solver.SharedMutationJoiningSolver import SharedMutationJoiningSolver

from run_512cell_baseline_comparison import (
    MISSING_STATE,
    balanced_split,
    build_signed_adjacency,
    evaluate_cassiopeia_tree,
    evaluate_networkx_tree,
    infer_sampled_triplets,
    local_search_maxcut,
    simulate_tree,
    solve_baseline,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a 512-cell sampled TMC benchmark with local SMJ refinement."
    )
    parser.add_argument("--depth", type=int, default=10)
    parser.add_argument("--sites", type=int, nargs="+", default=[120])
    parser.add_argument("--missing-rates", type=float, nargs="+", default=[0.1])
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--mutation-rate", type=float, default=0.8)
    parser.add_argument("--triplet-samples", type=int, default=300000)
    parser.add_argument("--accuracy-samples", type=int, default=20000)
    parser.add_argument("--tmc-restarts", type=int, default=8)
    parser.add_argument("--refine-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=5122026)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/cassiopeia_512cell_refined_tmc_baseline_comparison.csv"),
    )
    return parser.parse_args()


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
            attach_smj_refinement(
                graph,
                parent,
                ground_tree,
                leaves,
                subset,
                prefix=parent,
            )
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


def run_condition(args, k_sites, missing_rate, repeat):
    seed = args.seed + 100000 * repeat + 1000 * k_sites + int(100 * missing_rate)
    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    start = time.perf_counter()
    tree = simulate_tree(k_sites, missing_rate, args.depth, args.mutation_rate)
    simulation_seconds = time.perf_counter() - start
    leaves = list(tree.leaves)

    start = time.perf_counter()
    sampled_triplets = infer_sampled_triplets(tree, args.triplet_samples, rng)
    triplet_seconds = time.perf_counter() - start

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

    rows = []
    start = time.perf_counter()
    refined_tree = refined_sampled_tmc_tree(
        tree,
        leaves,
        sampled_triplets,
        rng,
        args.tmc_restarts,
        args.refine_size,
    )
    solve_seconds = time.perf_counter() - start
    eval_start = time.perf_counter()
    triplet_accuracy, rf_similarity = evaluate_networkx_tree(
        tree,
        refined_tree,
        args.accuracy_samples,
    )
    eval_seconds = time.perf_counter() - eval_start
    rows.append(
        base
        | {
            "algorithm": f"SampledRecursiveTMC_SMJ{args.refine_size}",
            "triplet_accuracy": triplet_accuracy,
            "rf_similarity": rf_similarity,
            "triplet_inference_seconds": triplet_seconds,
            "solve_seconds": solve_seconds,
            "evaluation_seconds": eval_seconds,
            "status": "ok",
            "error": "",
        }
    )
    print(
        f"done SampledRecursiveTMC_SMJ{args.refine_size} cells={len(leaves)} "
        f"sites={k_sites} missing={missing_rate:.2f} triplet_acc={triplet_accuracy:.3f} "
        f"rf={rf_similarity:.3f} sampled={len(sampled_triplets)} solve_sec={solve_seconds:.1f}",
        flush=True,
    )

    for algorithm in ("VanillaGreedy", "UPGMA", "NeighborJoining"):
        start = time.perf_counter()
        try:
            recon = solve_baseline(tree, algorithm)
            baseline_solve_seconds = time.perf_counter() - start
            eval_start = time.perf_counter()
            baseline_triplet, baseline_rf = evaluate_cassiopeia_tree(
                tree,
                recon,
                args.accuracy_samples,
            )
            baseline_eval = time.perf_counter() - eval_start
            status = "ok"
            error = ""
        except Exception as exc:
            baseline_solve_seconds = time.perf_counter() - start
            baseline_eval = 0.0
            baseline_triplet = np.nan
            baseline_rf = np.nan
            status = "error"
            error = f"{type(exc).__name__}: {exc}"
        rows.append(
            base
            | {
                "algorithm": algorithm,
                "triplet_accuracy": baseline_triplet,
                "rf_similarity": baseline_rf,
                "triplet_inference_seconds": 0.0,
                "solve_seconds": baseline_solve_seconds,
                "evaluation_seconds": baseline_eval,
                "status": status,
                "error": error,
            }
        )
        print(
            f"done {algorithm} cells={len(leaves)} sites={k_sites} missing={missing_rate:.2f} "
            f"triplet_acc={baseline_triplet:.3f} rf={baseline_rf:.3f} "
            f"solve_sec={baseline_solve_seconds:.1f} status={status}",
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
            f"{row['algorithm']:24s} cells={row['n_cells']} sites={row['n_sites']} "
            f"missing={row['missing_rate']:.2f} triplet={row['triplet_accuracy']:.3f} "
            f"rf={row['rf_similarity']:.3f} solve_sec={row['solve_seconds']:.1f} "
            f"status={row['status']}"
        )


if __name__ == "__main__":
    main()

