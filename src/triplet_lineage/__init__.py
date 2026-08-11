"""Utilities for triplet-based molecular lineage reconstruction."""

from importlib import import_module

_LAZY_EXPORTS = {
    "COLORS": ("plotting", "COLORS"),
    "configure_journal_style": ("plotting", "configure_journal_style"),
    "log10_formatter": ("plotting", "log10_formatter"),
    "compute_delta_star": ("theory", "compute_delta_star"),
    "sample_complexity_bound": ("theory", "sample_complexity_bound"),
    "solve_delta_star": ("theory", "solve_delta_star"),
    "calculate_rf_for_maxcut": ("metrics", "calculate_rf_for_maxcut"),
    "calculate_triplets_correct": ("metrics", "calculate_triplets_correct"),
    "find_triplet_structure": ("metrics", "find_triplet_structure"),
    "robinson_foulds_score": ("metrics", "robinson_foulds_score"),
    "robinson_foulds_zero": ("metrics", "robinson_foulds_zero"),
    "triplet_accuracy": ("metrics", "triplet_accuracy"),
    "build_tree_from_triplet_partition": (
        "reconstruction",
        "build_tree_from_triplet_partition",
    ),
    "build_tree_from_triplets": ("reconstruction", "build_tree_from_triplets"),
    "construct_triplet_graph": ("reconstruction", "construct_triplet_graph"),
    "find_recon_triplets": ("reconstruction", "find_recon_triplets"),
    "goemans_williamson_solve_desired_triplets": (
        "reconstruction",
        "goemans_williamson_solve_desired_triplets",
    ),
    "gw_partition": ("reconstruction", "gw_partition"),
    "hamming_distance": ("reconstruction", "hamming_distance"),
    "infer_triplets_from_mutations": ("reconstruction", "infer_triplets_from_mutations"),
    "percolation_solve": ("reconstruction", "percolation_solve"),
    "shared_mutation_solve": ("reconstruction", "shared_mutation_solve"),
    "complete_binary_missing_tree_sim": (
        "simulation",
        "complete_binary_missing_tree_sim",
    ),
    "complete_binary_topology_sim": ("simulation", "complete_binary_topology_sim"),
    "complete_binary_tree_sim": ("simulation", "complete_binary_tree_sim"),
    "depth_isomorphism": ("simulation", "depth_isomorphism"),
    "exponential_plus_c_topology_sim": (
        "simulation",
        "exponential_plus_c_topology_sim",
    ),
    "exponential_plus_c_tree_sim": ("simulation", "exponential_plus_c_tree_sim"),
    "overlay_mut_data": ("simulation", "overlay_mut_data"),
    "overlay_mutation_data": ("simulation", "overlay_mutation_data"),
}


def __getattr__(name):
    """Load Cassiopeia-dependent helpers only when they are requested."""
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module 'triplet_lineage' has no attribute {name!r}")

    module_name, attribute_name = _LAZY_EXPORTS[name]
    module = import_module(f".{module_name}", __name__)
    value = getattr(module, attribute_name)
    globals()[name] = value
    return value


__all__ = [
    "COLORS",
    "build_tree_from_triplet_partition",
    "build_tree_from_triplets",
    "calculate_rf_for_maxcut",
    "calculate_triplets_correct",
    "complete_binary_missing_tree_sim",
    "complete_binary_topology_sim",
    "complete_binary_tree_sim",
    "compute_delta_star",
    "configure_journal_style",
    "construct_triplet_graph",
    "depth_isomorphism",
    "exponential_plus_c_topology_sim",
    "exponential_plus_c_tree_sim",
    "find_recon_triplets",
    "find_triplet_structure",
    "goemans_williamson_solve_desired_triplets",
    "gw_partition",
    "hamming_distance",
    "infer_triplets_from_mutations",
    "log10_formatter",
    "overlay_mut_data",
    "overlay_mutation_data",
    "percolation_solve",
    "robinson_foulds_score",
    "robinson_foulds_zero",
    "sample_complexity_bound",
    "shared_mutation_solve",
    "solve_delta_star",
    "triplet_accuracy",
]
