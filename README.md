

# Code and Simulations for Triplet-Based Lineage Inference

This repository provides code and simulation tools for studying triplet-based lineage tree reconstruction from molecular recording data.  
This project uses the Cassiopeia lineage tracing framework, which is licensed under the MIT License by Yosef Lab. The analysis makes use of the Cassiopeia package for lineage tracing, available [here](https://github.com/YosefLab/Cassiopeia).

This codebase accompanies the manuscript:

**Theoretical Limits and Design Principles for Reliable Molecular Lineage Recorders**

The project focuses on theoretical analysis and simulation-based evaluation of triplet-based inference methods, including Max-Cut–based reconstruction algorithms, under realistic noise and sparsity constraints.

---

## Overview

Genetically encoded molecular recorders encode lineage history into heritable molecular barcodes that accumulate edits over time. In practice, these barcodes are sparse, noisy, and subject to biological and technical constraints, which limits reconstruction accuracy.

This repository provides:

- A theoretical framework relating local triplet error rates to global lineage tree reconstruction accuracy  
- Analytical bounds on reconstruction accuracy and minimum barcode length under simplified noise models  
- Simulation code for validating theoretical predictions under realistic recording noise  
- Implementations of Triplet Max-Cut (TMC) and baseline lineage reconstruction algorithms  

All code is intended to reproduce the simulations and figures reported in the accompanying manuscript.

---

## Repository Structure

```text
.
├── theory/
│   └── parameter_scaling.py              # Scaling laws for sample complexity and design parameters
│
├── results/
│   ├── rf_scores_with_missing.csv        # RF distances for baseline methods with missing data
│   ├── rf_scores_without_missing.csv     # RF distances for baseline methods without missing data
│   ├── triplet_scores_with_missing.csv   # Triplet accuracy under missing data
│   └── triplet_scores_without_missing.csv# Triplet accuracy without missing data
│
├── utilities/
│   ├── algorithm.py                     # Triplet-based reconstruction algorithms (TMC and variants)
│   ├── plot_style.py                    # Plotting utilities and shared figure style
│   ├── generate_tree.py                 # Lineage tree and barcode simulation utilities(from cassiopeia)
│   └── simulation_score.py              # Evaluation metrics for tree and triplet reconstruction
│
├── notebooks/
│   ├── validate_bounds.ipynb             # Empirical validation of theoretical accuracy bounds
│   ├── parameter_sensitivity.ipynb       # Sensitivity analysis over design parameters
│   ├── violin_plot.ipynb                 # Distributional analysis of reconstruction errors
│   └── plot_style.py                     # Notebook-specific plotting configuration
│
├── figures/
│   └── *.pdf                             # Figures generated from notebooks
│
├── requirements.txt
└── README.md

```

## Installation

### Requirements

- Python ≥ 3.8  
- NumPy  
- SciPy  
- NetworkX  
- Matplotlib  

Install all dependencies using:

```bash
pip install -r requirements.txt
```

## Usage

### 1. Simulate lineage trees and barcodes

### 2. Run Triplet Max-Cut reconstruction

### 3. Validate theoretical accuracy bounds

```
python notebooks/validate_bounds.ipynb
```

# License

This project is released under the MIT License.

---

#  Contact

Yuting Lin,
School of Mathematics,
Sun Yat-sen University,
Email: linyt66@mail2.sysu.edu.cn,
