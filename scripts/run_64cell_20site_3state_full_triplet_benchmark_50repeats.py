"""Run a 64-cell full-triplet benchmark with 5 edit states and 50 topology repeats."""

from __future__ import annotations

import csv
import itertools
import time
from pathlib import Path

import numpy as np

from run_512cell_balanced_stress_benchmark import (
    balanced_random_topology,
    evaluate_cassiopeia_tree,
    evaluate_networkx_tree,
    observed_hamming_rate,
    overlay_cas9_data,
    refined_tmc_tree,
    solve_baseline,
)


def infer_all_confident_triplets(tree, min_signal_gap: float):
    leaves = list(tree.leaves)
    values = tree.character_matrix.loc[leaves].to_numpy()
    inferred = []
    for a, b, c in itertools.combinations(range(len(leaves)), 3):
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


def run_repeat(repeat: int, missing_rate: float):
    seed = 6403005 + int(missing_rate * 1000) + repeat
    rng = np.random.default_rng(seed)
    np.random.seed(seed)
    tree = balanced_random_topology(64, rng)
    overlay_cas9_data(
        tree,
        n_sites=20,
        n_edit_states=3,
        dominant_state_prob=0.0,
        mutation_rate=0.8,
        missing_rate=missing_rate,
    )
    leaves = list(tree.leaves)
    rows = []

    start = time.perf_counter()
    triplets = infer_all_confident_triplets(tree, min_signal_gap=0.03)
    triplet_seconds = time.perf_counter() - start

    start = time.perf_counter()
    tmc_graph = refined_tmc_tree(tree, leaves, triplets, rng, restarts=4, refine_size=16)
    solve_seconds = time.perf_counter() - start
    triplet_accuracy, rf_similarity = evaluate_networkx_tree(tree, tmc_graph, accuracy_samples=5000)
    rows.append({
        "repeat": repeat,
        "seed": seed,
        "n_cells": 64,
        "n_sites": 20,
        "edit_states": 3,
        "missing_rate": missing_rate,
        "algorithm": "FullTripletTMC_SMJ16",
        "possible_triplets": 41664,
        "inferred_triplets": len(triplets),
        "triplet_coverage": len(triplets) / 41664,
        "triplet_accuracy": triplet_accuracy,
        "rf_similarity": rf_similarity,
        "triplet_inference_seconds": triplet_seconds,
        "solve_seconds": solve_seconds,
    })
    print(
        f"repeat={repeat} missing={missing_rate:.1f} TMC triplet={triplet_accuracy:.3f} "
        f"rf={rf_similarity:.3f} inferred={len(triplets)} infer={triplet_seconds:.2f}s solve={solve_seconds:.2f}s",
        flush=True,
    )

    for algorithm in ("VanillaGreedy", "UPGMA", "NeighborJoining"):
        start = time.perf_counter()
        recon = solve_baseline(tree, algorithm)
        solve_seconds = time.perf_counter() - start
        triplet_accuracy, rf_similarity = evaluate_cassiopeia_tree(tree, recon, accuracy_samples=5000)
        rows.append({
            "repeat": repeat,
            "seed": seed,
            "n_cells": 64,
            "n_sites": 20,
            "edit_states": 3,
            "missing_rate": missing_rate,
            "algorithm": algorithm,
            "possible_triplets": 41664,
            "inferred_triplets": len(triplets),
            "triplet_coverage": len(triplets) / 41664,
            "triplet_accuracy": triplet_accuracy,
            "rf_similarity": rf_similarity,
            "triplet_inference_seconds": 0.0,
            "solve_seconds": solve_seconds,
        })
        print(
            f"repeat={repeat} missing={missing_rate:.1f} {algorithm} triplet={triplet_accuracy:.3f} "
            f"rf={rf_similarity:.3f} solve={solve_seconds:.2f}s",
            flush=True,
        )
    return rows


def main() -> None:
    rows = []
    for missing_rate in (0.0, 0.3):
        for repeat in range(50):
            rows.extend(run_repeat(repeat, missing_rate))
    out = Path("results/cassiopeia_64cell_20site_3state_full_triplets_50repeats.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
