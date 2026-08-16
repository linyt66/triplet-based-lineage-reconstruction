"""Summarize the 64-cell full-triplet Cassiopeia benchmark.

The script aggregates replicate-level results by missing-data condition and
algorithm. It is intentionally small so it can be reused after extending the
number of repeats for publication figures.
"""

from pathlib import Path

import pandas as pd


RESULTS_DIR = Path("results")
INPUT_CSV = RESULTS_DIR / "cassiopeia_64cell_30site_3state_full_triplets_50repeats.csv"
OUTPUT_CSV = RESULTS_DIR / "cassiopeia_64cell_30site_3state_full_triplets_50repeats_summary.csv"


def main() -> None:
    """Write a compact per-condition summary table."""
    df = pd.read_csv(INPUT_CSV)
    summary = (
        df.groupby(["missing_rate", "algorithm"])
        .agg(
            n=("repeat", "count"),
            triplet_mean=("triplet_accuracy", "mean"),
            triplet_median=("triplet_accuracy", "median"),
            rf_mean=("rf_similarity", "mean"),
            rf_median=("rf_similarity", "median"),
            infer_mean=("triplet_inference_seconds", "mean"),
            solve_mean=("solve_seconds", "mean"),
            inferred_triplets_mean=("inferred_triplets", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(OUTPUT_CSV, index=False)
    print(summary.round(4).to_string(index=False))
    print(f"Wrote {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
