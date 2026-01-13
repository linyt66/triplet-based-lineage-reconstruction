

# Code and Simulations for Triplet-Based Lineage Inference

This repository provides code and simulation tools for studying triplet-based lineage tree reconstruction from molecular recording data.  
The analysis makes use of the Cassiopeia package for lineage tracing, available [here](https://github.com/YosefLab/Cassiopeia).

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
│   ├── triplet_error_bounds.py      # Numerical evaluation of theoretical bounds
│   └── parameter_scaling.py         # Scaling behavior of design parameters
│
├── simulation/
│   ├── generate_tree.py             # Lineage tree simulation with barcode evolution
│
├── results/
│   └── *.csv                        # UPGMA, NJ, and greedy evaluation
│ 
├── algorithms/
│   ├── tmc.py                       # Triplet Max-Cut (TMC) implementation
│   ├── oto.py                       # One-Time Optimization (OTO)
│   └── baselines.py                 # UPGMA, NJ, and greedy heuristics
│
├── notebooks/
│   ├── validate_bounds.ipynb        # Comparison of theoretical and empirical error rates
│   ├── parameter_sensitivity.ipynb  # Scaling with λ, ℓ_min, p_miss, p_share
│   └── rf_distance_eval.ipynb       # Robinson–Foulds distance evaluation
│
├── figures/
│   └── *.pdf                        # Generated figures
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
