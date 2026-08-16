# Triplet-Based Lineage Reconstruction

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Paper](https://img.shields.io/badge/Paper-Pending_Publication-green.svg)](#)

This repository contains the code, simulations, precomputed results, and
publication figures for the manuscript:

> **Theoretical Limits and Experimental Design Principles for Single-Cell
> Molecular Lineage Recording**

The project studies triplet-based inference for molecular lineage recording data.
It introduces the Optimal Triplet Oracle (OTO) framework and implements Triplet
Max-Cut (TMC), a recursive graph-partitioning reconstruction algorithm for noisy,
sparse, and homoplasy-prone molecular barcodes.

The code uses the Cassiopeia lineage tracing framework from Yosef Lab:
<https://github.com/YosefLab/Cassiopeia>.

## Repository Layout

```text
.
|-- LICENSE
|-- README.md
|-- pyproject.toml
|-- requirements.txt
|-- src/
|   `-- triplet_lineage/
|       |-- __init__.py
|       |-- metrics.py
|       |-- plotting.py
|       |-- reconstruction.py
|       |-- simulation.py
|       `-- theory.py
|-- notebooks/
|   |-- algorithm_performance.ipynb
|   |-- algorithm_performancev2.ipynb
|   |-- figure-6.ipynb
|   |-- validate_bounds.ipynb
|   `-- violin_plot.ipynb
|-- results/
|   |-- rf_scores_with_missing.csv
|   |-- rf_scores_without_missing.csv
|   |-- triplet_scores_with_missing.csv
|   `-- triplet_scores_without_missing.csv
`-- figures/
    |-- capacity-scaling-symlog.pdf
    |-- figure-01.pdf
    |-- figure-02.pdf
    |-- figure-03.pdf
    |-- figure-04.pdf
    |-- figure-05.pdf
    `-- reliability-bound.pdf
```

## Installation

Create an isolated Python environment, then install the repository in editable
mode:

```bash
git clone https://github.com/linyt66/triplet-based-lineage-reconstruction.git
cd triplet-based-lineage-reconstruction
pip install -e .
```

Alternatively, install the pinned runtime dependencies directly:

```bash
pip install -r requirements.txt
```

For development and validation, install the optional test dependency:

```bash
pip install -e ".[dev]"
```

## Code Modules

- `triplet_lineage.simulation`: complete binary and asynchronous topology
  simulators, Cas9 mutation overlays, and missing-data simulations.
- `triplet_lineage.reconstruction`: triplet extraction, signed graph
  construction, exact small-graph Max-Cut, recursive TMC reconstruction, Shared
  Mutation Joining, and percolation baselines.
- `triplet_lineage.metrics`: triplet accuracy and Robinson-Foulds evaluation.
- `triplet_lineage.theory`: OTO/TMC accuracy bounds, admissible triplet-error
  inversion, binomial reliability calculations, and site-complexity helpers.
- `triplet_lineage.plotting`: journal-style matplotlib defaults and shared color
  palette.

## Reproducing Experiments

Launch Jupyter from the repository root:

```bash
jupyter notebook
```

Recommended reproduction order:

1. Run `notebooks/algorithm_performance.ipynb` to simulate lineage trees,
   reconstruct them with TMC and baseline methods, and regenerate the CSV files
   in `results/`.
2. Run `notebooks/validate_bounds.ipynb` to compare empirical accuracy with the
   OTO and ITO theoretical limits.
3. Run `notebooks/algorithm_performancev2.ipynb` to estimate empirical recording
   capacity and reliability curves.
4. Run `notebooks/violin_plot.ipynb` and `notebooks/figure-6.ipynb` to
   regenerate the publication figures in `figures/`.

The precomputed CSV files in `results/` are included so that plotting notebooks
can be run without repeating the full simulation sweep.

## Suggested Validation Extensions

The unit tests intentionally remain lightweight. For manuscript validation,
prefer notebook or script-level experiments that can run for longer and save
summary CSV files:

- **Recursive TMC ablation:** compare recursive TMC against a single-level
  Max-Cut split followed by Shared Mutation Joining to isolate the effect of the
  recursive top-down step used in the manuscript.
- **Dropout stress test:** sweep `p_miss` from 0.0 to 0.3 and report both triplet
  accuracy and normalized RF similarity, matching the manuscript reliability
  claims.
- **Triplet abstention analysis:** vary the `min_signal_gap` threshold in
  `infer_triplets_from_mutations` and report coverage versus accuracy for
  resolvable triplets.
- **Theory-to-simulation calibration:** use `triplet_lineage.theory` to compute
  OTO/TMC lower bounds and overlay them with empirical mean and confidence
  intervals from repeated simulations.

## Testing

Run the lightweight test suite from the repository root:

```bash
python -m unittest discover -s tests -v
```

The tests cover package import behavior, theoretical helper functions, plotting
utilities, repository hygiene, triplet graph construction, dropout-aware triplet
inference, and recursive TMC reconstruction. Tests that require optional
scientific dependencies are skipped automatically when those dependencies are not
installed.

If `pytest` is installed, the same tests can also be run with:

```bash
pytest
```

## License

This project is released under the MIT License. See `LICENSE` for details.

Cassiopeia is distributed under its own MIT License by Yosef Lab.
