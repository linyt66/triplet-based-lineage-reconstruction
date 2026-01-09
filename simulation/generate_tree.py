
"""This file stores the utilities used in the simulations of the Threshold and 
Bottom-Up Algorithms. Requires the Cassiopeia software package: 
https://github.com/YosefLab/Cassiopeia."""

from typing import Dict, Tuple
import networkx as nx
import copy
import random
from itertools import combinations
from scipy.linalg import eigh
import cvxgraphalgs as cvxgr
from cvxgraphalgs.structures import Cut
import inspect
import networkx as nx
import pandas as pd
import numpy as np
import pickle as pic
import copy
import random
import itertools
import cassiopeia as cass
import cassiopeia.critique as critique
import cassiopeia.data as data
import cassiopeia.simulator as simulator
from cassiopeia.simulator.TreeSimulator import TreeSimulatorError
from cassiopeia.solver.VanillaGreedySolver import VanillaGreedySolver
from cassiopeia.solver.PercolationSolver import PercolationSolver
from cassiopeia.solver.SharedMutationJoiningSolver import SharedMutationJoiningSolver
from cassiopeia.solver import dissimilarity_functions

def complete_binary_tree_sim(k_cand: int, q_dist: Dict[int, float], lamb: float, depth: int) -> data.CassiopeiaTree:
    """Simulates a lineage tracing experiment with a complete binary topology.

    Args:
        k_cand: The number of characters
        q_dist: The state distribution
        lamb: The value of lambda
        depth: The depth of the tree (including the unary implicit root edge)

    Returns:
        A CassiopeiaTree containing the lineage tracing experiment
    """

    # Creating the binary networkx tree
    tree = nx.balanced_tree(2, (depth - 1), create_using = nx.DiGraph)
    tree = nx.relabel_nodes(tree, dict(zip(list(tree.nodes), [i + 1 for i in tree.nodes])))
     
    # Appending the implicit root edge
    tree.add_edge(0, 1)
    # Storing the networkx tree in a CassiopeiaTree container object
    complete_binary = data.CassiopeiaTree(tree = tree)

    # Relabling the leaves
    leaves = complete_binary.leaves
    complete_binary.relabel_nodes(dict(zip(leaves, ['c' + str(i) for i in leaves])))
    
    # Normalizing the times and edge lengths by the total time of the tree
    time_dict = {}
    for i in complete_binary.nodes:
        time_dict[i] = complete_binary.get_time(i)/depth
    complete_binary.set_times(time_dict)
    
    # Generating the lineage tracing data with the supplied parameters and
    # with no missing data
    lt_sim = simulator.Cas9LineageTracingDataSimulator(
        number_of_cassettes = k_cand,
        size_of_cassette = 1,
        mutation_rate = lamb,
        state_priors = q_dist,
        heritable_silencing_rate = 0,
        stochastic_silencing_rate = 0,
        heritable_missing_data_state = -1,
        stochastic_missing_data_state = -1,
    )
    lt_sim.overlay_data(complete_binary)

    return complete_binary

def exponential_plus_c_tree_sim(k_cand: int, q_dist: Dict[int, float], lamb: float, depth: int) -> data.CassiopeiaTree:
    """Simulates a lineage tracing experiment with an asynchronous topology.

    Args:
        k_cand: The number of characters
        q_dist: The state distribution
        lamb: The value of lambda
        depth: The depth of the tree (including the unary implicit root edge)

    Returns:
        A CassiopeiaTree containing the lineage tracing experiment
    """
    # We define the number of leaves wanted as 2**(depth - 1)
    num_cells = 2**(depth - 1)
    size = 0
    c = 0.05

    # The birth rate and death rate are hard-coded such that they produce
    # trees with an average of 256 leaves
    birth_rate = 23.697339506322916 
    death_rate = 2.3697339506322916

    # The waiting distributions for our Bellman-Harris process with death
    initial_birth_scale = 1/birth_rate
    birth_waiting_distribution = lambda scale: np.random.exponential(scale) + c
    death_waiting_distribution = lambda: np.random.exponential(1/death_rate + c)

    # We generate trees until the number of leaves is between a factor of 0.8
    # and 1.2 of our wanted number of leaves
    while size < num_cells * 0.8 or size > num_cells * 1.2:
        try:
            bd_sim = simulator.BirthDeathFitnessSimulator(
                birth_waiting_distribution = birth_waiting_distribution,
                initial_birth_scale = initial_birth_scale,
                death_waiting_distribution = death_waiting_distribution,
                experiment_time = 1,
            )
            topology = bd_sim.simulate_tree()
            size = topology.n_cell
        except TreeSimulatorError:
            size = 0

    # Relabeling the leaves
    topology.relabel_nodes(dict(zip(topology.leaves, ["c" + i for i in topology.leaves])))

    # Generating the lineage tracing data with the supplied parameters and
    # with no missing data
    lt_sim = simulator.Cas9LineageTracingDataSimulator(
        number_of_cassettes = k_cand,
        size_of_cassette = 1,
        mutation_rate = lamb,
        state_priors = q_dist,
        heritable_silencing_rate = 0,
        stochastic_silencing_rate = 0,
        heritable_missing_data_state = -1,
        stochastic_missing_data_state = -1,
    )
    lt_sim.overlay_data(topology)

    return topology

def complete_binary_missing_tree_sim(k_cand: int, q_dist: Dict[int, float], lamb: float, depth: int) -> data.CassiopeiaTree:
    """Simulates a lineage tracing experiment with a complete binary topology
    and missing data.

    Args:
        k_cand: The number of characters
        q_dist: The state distribution
        lamb: The value of lambda
        depth: The depth of the tree (including the unary implicit root edge)

    Returns:
        A CassiopeiaTree containing the lineage tracing experiment
    """

    # Creating the binary networkx tree
    tree = nx.balanced_tree(2, (depth - 1), create_using = nx.DiGraph)
    tree = nx.relabel_nodes(tree, dict(zip(list(tree.nodes), [i + 1 for i in tree.nodes])))

    # Appending the implicit root edge
    tree.add_edge(0, 1)

    # Storing the networkx tree in a CassiopeiaTree container object
    complete_binary = data.CassiopeiaTree(tree = tree)

    # Relabling the leaves
    leaves = complete_binary.leaves
    complete_binary.relabel_nodes(dict(zip(leaves, ['c' + str(i) for i in leaves])))

    # Normalizing the times and edge lengths by the total time of the tree
    time_dict = {}
    for i in complete_binary.nodes:
        time_dict[i] = complete_binary.get_time(i)/depth
    complete_binary.set_times(time_dict)
    
    # Generating the lineage tracing data with the supplied parameters and
    # with 10% stochastic missing data
    lt_sim = simulator.Cas9LineageTracingDataSimulator(
        number_of_cassettes = k_cand,
        size_of_cassette = 1,
        mutation_rate = lamb,
        state_priors = q_dist,
        heritable_silencing_rate = 0,
        stochastic_silencing_rate = 0.1,
        heritable_missing_data_state = -1,
        stochastic_missing_data_state = -1,
    )
    lt_sim.overlay_data(complete_binary)

    return complete_binary

def complete_binary_missing_tree_sim(k_cand: int, q_dist: Dict[int, float], lamb: float, depth: int) -> data.CassiopeiaTree:
    """Simulates a lineage tracing experiment with a complete binary topology
    and missing data.

    Args:
        k_cand: The number of characters
        q_dist: The state distribution
        lamb: The value of lambda
        depth: The depth of the tree (including the unary implicit root edge)

    Returns:
        A CassiopeiaTree containing the lineage tracing experiment
    """

    # Creating the binary networkx tree
    tree = nx.balanced_tree(2, (depth - 1), create_using = nx.DiGraph)
    tree = nx.relabel_nodes(tree, dict(zip(list(tree.nodes), [i + 1 for i in tree.nodes])))

    # Appending the implicit root edge
    tree.add_edge(0, 1)

    # Storing the networkx tree in a CassiopeiaTree container object
    complete_binary = data.CassiopeiaTree(tree = tree)

    # Relabling the leaves
    leaves = complete_binary.leaves
    complete_binary.relabel_nodes(dict(zip(leaves, ['c' + str(i) for i in leaves])))

    # Normalizing the times and edge lengths by the total time of the tree
    time_dict = {}
    for i in complete_binary.nodes:
        time_dict[i] = complete_binary.get_time(i)/depth
    complete_binary.set_times(time_dict)
    
    # Generating the lineage tracing data with the supplied parameters and
    # with 10% stochastic missing data
    lt_sim = simulator.Cas9LineageTracingDataSimulator(
        number_of_cassettes = k_cand,
        size_of_cassette = 1,
        mutation_rate = lamb,
        state_priors = q_dist,
        heritable_silencing_rate = 0,
        stochastic_silencing_rate = 0.1,
        heritable_missing_data_state = -1,
        stochastic_missing_data_state = -1,
    )
    lt_sim.overlay_data(complete_binary)

    return complete_binary
