"""Evaluation metrics for triplet-based lineage reconstruction."""

from typing import Tuple

import cassiopeia.critique as critique
import cassiopeia.data as data
import networkx as nx
import numpy as np


def find_triplet_structure(triplet: Tuple[str, str, str], tree: nx.DiGraph) -> str:
    """Return the resolved ingroup of a rooted triplet, or '-' for a tie."""
    a, b, c = triplet

    ancestors = {
        a: set(nx.ancestors(tree, a)),
        b: set(nx.ancestors(tree, b)),
        c: set(nx.ancestors(tree, c)),
    }
    shared_ancestor_counts = {
        "ab": len(ancestors[a] & ancestors[b]),
        "ac": len(ancestors[a] & ancestors[c]),
        "bc": len(ancestors[b] & ancestors[c]),
    }

    best = max(shared_ancestor_counts, key=shared_ancestor_counts.get)
    values = list(shared_ancestor_counts.values())
    return best if values.count(shared_ancestor_counts[best]) == 1 else "-"


def calculate_triplets_correct(
    ground_tree: data.CassiopeiaTree,
    recon_tree: nx.DiGraph,
    n_samples: int = 5000,
) -> float:
    """Estimate the fraction of randomly sampled triplets resolved correctly."""
    ground_topology = ground_tree.get_tree_topology()
    leaves = list(ground_tree.leaves)

    correct = 0
    for _ in range(n_samples):
        sampled_triplet = tuple(np.random.choice(leaves, 3, replace=False))
        correct += int(
            find_triplet_structure(sampled_triplet, recon_tree)
            == find_triplet_structure(sampled_triplet, ground_topology)
        )

    return correct / n_samples if n_samples else 0.0


def triplet_accuracy(
    ground_tree: data.CassiopeiaTree,
    recon_tree: nx.DiGraph,
    n_samples: int = 5000,
) -> float:
    """Alias for ``calculate_triplets_correct``."""
    return calculate_triplets_correct(ground_tree, recon_tree, n_samples)


def robinson_foulds_score(
    ground_tree: data.CassiopeiaTree,
    recon_tree: data.CassiopeiaTree,
) -> bool:
    """Return True when the Robinson-Foulds distance is zero."""
    rf_distance, _ = critique.compare.robinson_foulds(ground_tree, recon_tree)
    return rf_distance == 0


def robinson_foulds_zero(
    ground_tree: data.CassiopeiaTree,
    recon_tree: data.CassiopeiaTree,
) -> bool:
    """Compatibility alias for ``robinson_foulds_score``."""
    return robinson_foulds_score(ground_tree, recon_tree)


def calculate_rf_for_maxcut(
    ground_tree: data.CassiopeiaTree,
    recon_tree: nx.DiGraph,
) -> Tuple[float, float]:
    """Compute raw and normalized RF distance for a NetworkX TMC tree."""
    cassiopeia_recon = data.CassiopeiaTree(tree=recon_tree)
    rf_distance, rf_max = critique.compare.robinson_foulds(
        ground_tree,
        cassiopeia_recon,
    )
    normalized = rf_distance / rf_max if rf_max else 0.0
    return rf_distance, normalized
