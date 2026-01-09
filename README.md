# Lineage Triplet Theory

This repository contains the code and simulation framework accompanying the manuscript:

**Fundamental Limits on Phylogeny Inference from Single-Cell Lineage Tracing**  
*(submitted to Systematic Biology)*

This project studies the theoretical limits and optimal experimental design for reconstructing lineage trees from genetically encoded molecular recording data, with a focus on triplet-based inference and Max-Cut–based reconstruction algorithms.

---

## Overview

Genetically encoded molecular recorders convert lineage history into heritable molecular barcodes that accumulate changes over time. However, these recordings are sparse, noisy, and constrained by biological and technical limitations. This repository provides:

- A theoretical framework linking local triplet error rates to global tree reconstruction accuracy  
- Provable bounds on reconstruction accuracy and minimum barcode length  
- Simulation code validating theoretical predictions under realistic noise models  
- Implementations of Triplet Max-Cut (TMC) and baseline reconstruction algorithms  

All code is designed to reproduce the simulations and figures reported in the manuscript.

---

## Repository Structure

```text
.
├── theory/
│   ├── triplet_error_bounds.py      # Numerical evaluation of theoretical bounds
│   └── parameter_scaling.py         # Scaling laws for design parameters
│
├── simulation/
│   ├── generate_tree.py             # Random lineage tree generation
│   ├── simulate_barcodes.py         # CRISPR-style barcode evolution
│   └── noise_models.py              # Missing data and shared-edit models
│
├── algorithms/
│   ├── tmc.py                       # Triplet Max-Cut (TMC) implementation
│   ├── oto.py                       # One-Time Optimization (OTO)
│   └── baselines.py                 # UPGMA, NJ, greedy heuristics
│
├── experiments/
│   ├── validate_bounds.py           # Reproduces Fig. 5 (theoretical vs empirical)
│   ├── parameter_sensitivity.py     # Scaling with λ, ℓ_min, p_miss, p_share
│   └── rf_distance_eval.py          # Robinson–Foulds distance analysis
│
├── figures/
│   └── *.pdf                        # Generated figures
│
├── requirements.txt
└── README.md

