# Code and Simulations for Triplet-Based Lineage Inference

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Paper](https://img.shields.io/badge/Paper-Pending_Publication-green.svg)](#)

This repository provides code and simulation tools for studying triplet-based lineage tree reconstruction from molecular recording data.

This project uses the **Cassiopeia** lineage tracing framework, which is licensed under the MIT License by Yosef Lab. The analysis heavily integrates the Cassiopeia package for lineage tracing simulation and data handling, available [here](https://github.com/YosefLab/Cassiopeia).

This codebase accompanies the manuscript:

> **Theoretical Limits and Experimental Design Principles for Single-Cell Molecular Lineage Recording**  
> *(Authors' Names, Journal/Conference, Year)*

The project focuses on the theoretical analysis and simulation-based evaluation of triplet-based inference methods—specifically, the **Triplet Max-Cut (TMC)** reconstruction algorithm—under realistic noise, homoplasy, and sparsity constraints.

---

# 📖 Overview

Genetically encoded molecular recorders capture lineage history into heritable molecular barcodes that accumulate edits over time. In practice, these barcodes are sparse, noisy, and subject to biological and technical constraints (e.g., mutational homoplasy and dropout), which limit global reconstruction accuracy.

This repository provides:

- A theoretical framework relating local triplet error rates to global lineage tree reconstruction accuracy.
- Analytical bounds on reconstruction accuracy and minimal barcode length under various noise models.
- Simulation code for validating theoretical predictions under realistic recording noise regimes.
- Implementations of the **Triplet Max-Cut (TMC)** algorithm alongside baseline lineage reconstruction heuristics (e.g., Neighbor-Joining and UPGMA).

All code is intended to fully reproduce the simulations, theoretical evaluations, and figures reported in the accompanying manuscript.

---

# 🗂 Repository Structure

```text
.
├── README.md
├── LICENSE.txt
├── requirements.txt
├── theory/
│   └── parameter_scaling.py
├── results/
│   ├── rf_scores_with_missing.csv
│   ├── rf_scores_without_missing.csv
│   ├── triplet_scores_with_missing.csv
│   └── triplet_scores_without_missing.csv
├── utilities/
│   ├── algorithm.py
│   ├── generate_tree.py
│   ├── simulation_score.py
│   ├── plot_style.py
│   └── utilities.py
├── notebooks/
│   ├── algorithm_performance.ipynb
│   ├── algorithm_performancev2.ipynb
│   ├── parameter_sensitivity.ipynb
│   ├── validate_bounds.ipynb
│   ├── violin_plot.ipynb
│   ├── figure-6.ipynb
│   └── plot_style.py
├── figures/
│   ├── Fig_Capacity_Final_Symlog.pdf
│   ├── Fig_Reliability_Final.pdf
│   ├── figure-1.pdf
│   ├── figure-2.pdf
│   ├── figure-3.pdf
│   ├── figure-4.pdf
│   ├── figure-5.pdf
│   └── figure-6.pdf
```

---

# ⚙️ Installation

We recommend using a virtual environment (e.g., `conda` or `venv`) to run this project.

## Requirements

- Python ≥ 3.8
- NumPy
- SciPy
- NetworkX
- Matplotlib
- Pandas

## Setup

Clone the repository and install the required dependencies:

```bash
git clone https://github.com/linyt66/triplet-based-lineage-reconstruction.git
cd triplet-based-lineage-reconstruction
pip install -r requirements.txt
```

---

# 🚀 Usage

The core experiments are organized into intuitive Jupyter notebooks.

Launch Jupyter from the project root:

```bash
jupyter notebook
```

## 1. Simulate Lineage Trees and Molecular Barcodes

The simulation framework is implemented in:

```text
utilities/generate_tree.py
```

To reproduce the generation of ground-truth lineage trees and noisy molecular recording matrices (including CRISPR edits and missing data), execute the initial sections of

```text
notebooks/algorithm_performance.ipynb
```

---

## 2. Run Triplet Max-Cut (TMC) Reconstruction

To benchmark TMC against standard phylogenetic heuristics (Neighbor-Joining and UPGMA):

1. Open `notebooks/algorithm_performance.ipynb`.
2. Execute all cells.
3. The notebook will:
   - simulate lineage trees,
   - generate molecular recording data,
   - reconstruct trees using TMC,
   - compare against NJ and UPGMA,
   - compute Triplet Scores and Robinson–Foulds distances.

The resulting metrics will be saved in

```text
results/
```

---

## 3. Validate Theoretical Accuracy Bounds

To compare empirical reconstruction performance with the proposed **Optimal Triplet Oracle (OTO)** and **Ideal Triplet Oracle (ITO)** theoretical limits:

Open

```text
notebooks/validate_bounds.ipynb
```

and execute all cells.

The notebook either

- loads the precomputed simulation results from `results/`, or
- regenerates simulations,

to reproduce the theoretical accuracy curves and barcode capacity scaling laws reported in the manuscript.

---

## 4. Reproduce Publication Figures

Publication-quality figures can be regenerated using the notebooks

- `violin_plot.ipynb`
- `parameter_sensitivity.ipynb`
- `figure-6.ipynb`

All plotting styles are centralized in

```text
utilities/plot_style.py
```

The generated figures are saved into

```text
figures/
```

# 📄 License

This project is released under the **MIT License**.

See:

```text
LICENSE.txt
```

This project also interfaces with the **Cassiopeia** lineage tracing framework, which is distributed under its own MIT License.

For details, visit:

https://github.com/YosefLab/Cassiopeia