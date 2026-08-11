"""Simulation helpers for molecular lineage reconstruction experiments."""

from typing import Dict

import copy

import cassiopeia.data as data
import cassiopeia.simulator as simulator
import networkx as nx
import numpy as np
from cassiopeia.simulator.TreeSimulator import TreeSimulatorError


def _normalize_node_times(tree: data.CassiopeiaTree, depth: int) -> None:
    """Scale Cassiopeia node times by the requested tree depth in place."""
    time_dict = {
        node: tree.get_time(node) / depth
        for node in tree.nodes
    }
    tree.set_times(time_dict)


def _cas9_simulator(
    k_cand: int,
    q_dist: Dict[int, float],
    lamb: float,
    stochastic_silencing_rate: float = 0.0,
) -> simulator.Cas9LineageTracingDataSimulator:
    """Create the Cas9 simulator used across the manuscript experiments."""
    return simulator.Cas9LineageTracingDataSimulator(
        number_of_cassettes=k_cand,
        size_of_cassette=1,
        mutation_rate=lamb,
        state_priors=q_dist,
        heritable_silencing_rate=0,
        stochastic_silencing_rate=stochastic_silencing_rate,
        heritable_missing_data_state=-1,
        stochastic_missing_data_state=-1,
    )


def complete_binary_topology_sim(depth: int) -> data.CassiopeiaTree:
    """Generate a complete binary topology with normalized node times."""
    topology = nx.balanced_tree(2, depth - 1, create_using=nx.DiGraph)
    topology = nx.relabel_nodes(
        topology,
        dict(zip(list(topology.nodes), [node + 1 for node in topology.nodes])),
    )
    topology.add_edge(0, 1)

    complete_binary = data.CassiopeiaTree(tree=topology)
    leaves = complete_binary.leaves
    complete_binary.relabel_nodes({
        leaf: "c" + str(leaf)
        for leaf in leaves
    })
    _normalize_node_times(complete_binary, depth)

    return complete_binary


def complete_binary_tree_sim(
    k_cand: int,
    q_dist: Dict[int, float],
    lamb: float,
    depth: int,
) -> data.CassiopeiaTree:
    """Simulate a complete binary lineage tracing experiment."""
    tree = complete_binary_topology_sim(depth)
    _cas9_simulator(k_cand, q_dist, lamb).overlay_data(tree)
    return tree


def complete_binary_missing_tree_sim(
    k_cand: int,
    q_dist: Dict[int, float],
    lamb: float,
    depth: int,
    missing_rate: float = 0.1,
) -> data.CassiopeiaTree:
    """Simulate a complete binary experiment with stochastic missing data."""
    tree = complete_binary_topology_sim(depth)
    _cas9_simulator(
        k_cand,
        q_dist,
        lamb,
        stochastic_silencing_rate=missing_rate,
    ).overlay_data(tree)
    return tree


def exponential_plus_c_topology_sim(depth: int) -> data.CassiopeiaTree:
    """Generate an asynchronous Bellman-Harris birth-death topology."""
    target_leaves = 2 ** (depth - 1)
    size = 0
    offset = 0.05

    birth_rate = 23.697339506322916
    death_rate = 2.3697339506322916
    initial_birth_scale = 1 / birth_rate

    def birth_waiting_distribution(scale):
        return np.random.exponential(scale) + offset

    def death_waiting_distribution():
        return np.random.exponential(1 / death_rate + offset)

    # Resample until the topology has approximately the intended leaf count.
    while size < target_leaves * 0.8 or size > target_leaves * 1.2:
        try:
            bd_sim = simulator.BirthDeathFitnessSimulator(
                birth_waiting_distribution=birth_waiting_distribution,
                initial_birth_scale=initial_birth_scale,
                death_waiting_distribution=death_waiting_distribution,
                experiment_time=1,
            )
            topology = bd_sim.simulate_tree()
            size = topology.n_cell
        except TreeSimulatorError:
            size = 0

    topology.relabel_nodes({
        leaf: "c" + str(leaf)
        for leaf in topology.leaves
    })

    return topology


def exponential_plus_c_tree_sim(
    k_cand: int,
    q_dist: Dict[int, float],
    lamb: float,
    depth: int,
) -> data.CassiopeiaTree:
    """Simulate an asynchronous lineage tracing experiment."""
    topology = exponential_plus_c_topology_sim(depth)
    _cas9_simulator(k_cand, q_dist, lamb).overlay_data(topology)
    return topology


def overlay_mutation_data(
    topology: data.CassiopeiaTree,
    k_cand: int,
    q_dist: Dict[int, float],
    lamb: float,
) -> data.CassiopeiaTree:
    """Overlay Cas9 mutation data onto an existing topology."""
    _cas9_simulator(k_cand, q_dist, lamb).overlay_data(topology)
    return topology


def depth_isomorphism(
    d: float,
    ground_tree: data.CassiopeiaTree,
    recon_tree: data.CassiopeiaTree,
) -> bool:
    """Check whether two reconstructed topologies agree up to depth ``d``."""
    tree_g = copy.copy(ground_tree)
    tree_r = copy.copy(recon_tree)

    tree_g.collapse_unifurcations()
    tree_r.collapse_unifurcations()

    def check_split(node_g, node_r, correct_flag):
        groups_g = {
            child: set(tree_g.leaves_in_subtree(child))
            for child in tree_g.children(node_g)
        }
        groups_r = {
            child: set(tree_r.leaves_in_subtree(child))
            for child in tree_r.children(node_r)
        }

        matched_nodes = [
            (left, right)
            for left, leaves_g in groups_g.items()
            for right, leaves_r in groups_r.items()
            if leaves_g == leaves_r
        ]

        if not matched_nodes:
            correct_flag[0] = False
            return

        for left, right in matched_nodes:
            is_internal = len(tree_g.leaves_in_subtree(left)) > 1
            if tree_g.get_attribute(left, "time") <= d and correct_flag[0] and is_internal:
                check_split(left, right, correct_flag)

    correct_flag = [True]
    check_split(tree_g.root, tree_r.root, correct_flag)
    return correct_flag[0]


overlay_mut_data = overlay_mutation_data
