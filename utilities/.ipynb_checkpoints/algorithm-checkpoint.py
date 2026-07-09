from typing import Dict, Tuple
import networkx as nx
import pandas as pd
import numpy as np
import pickle as pic
import copy
import random
from itertools import combinations
from scipy.linalg import eigh
import itertools
import utilities
import cassiopeia.critique as critique
import cassiopeia.data as data
import cassiopeia.simulator as simulator
from cassiopeia.simulator.TreeSimulator import TreeSimulatorError
from cassiopeia.solver.VanillaGreedySolver import VanillaGreedySolver
from cassiopeia.solver.PercolationSolver import PercolationSolver
from cassiopeia.solver.SharedMutationJoiningSolver import SharedMutationJoiningSolver
from cassiopeia.solver import dissimilarity_functions
import cvxgraphalgs as cvxgr
from cvxgraphalgs.structures import Cut
import inspect


def construct_graph(triplets):
    """
    Construct a weighted graph based on triplet constraints.

    Args:
        triplets: List of triplets, where each triplet is a tuple of three node labels.

    Returns:
        G_simple: A weighted simple graph constructed from the triplets.
    """

    G = nx.MultiGraph()

    for triplet in triplets:
        a, b, c = triplet

        # Add edges according to TD-triplet constraints
        G.add_edge(a, b, weight=-2)
        G.add_edge(a, c, weight=1)
        G.add_edge(b, c, weight=1)

    # Convert the multigraph into a simple graph and accumulate edge weights
    G_simple = nx.Graph()

    for u, v, data in G.edges(data=True):
        if G_simple.has_edge(u, v):
            G_simple[u][v]['weight'] += data['weight']
        else:
            G_simple.add_edge(u, v, weight=data['weight'])

    return G_simple


def goemans_williamson_solve_desired_triplets(triplets):
    """
    Apply the Goemans-Williamson SDP relaxation to compute a bipartition
    that maximizes the satisfied triplet constraints.
    """

    G = construct_graph(triplets)

    sdp_cut = cvxgr.algorithms.goemans_williamson_weighted(G)

    best_partition = [set(), set()]
    best_partition[0] = sdp_cut.left
    best_partition[1] = sdp_cut.right

    S = set(best_partition[0])
    S_bar = set(best_partition[1])

    obeyed_triplets = []
    disobeyed_triplets = []

    for triplet in triplets:
        a, b, c = triplet

        if (a in S and b in S and c not in S) or (a not in S and b not in S and c in S):
            obeyed_triplets.append(triplet)
        else:
            disobeyed_triplets.append(triplet)

    cut_value = sdp_cut.evaluate_cut_size(G)

    return S, S_bar


def build_tree_from_triplet_partition(tree: data.CassiopeiaTree, triplets) -> nx.DiGraph:
    """
    Build a reconstructed tree from a triplet-based partition.

    Args:
        tree: CassiopeiaTree containing the character matrix
        triplets: Triplet constraints

    Returns:
        recon_tree: A reconstructed directed tree
    """

    S, S_bar = goemans_williamson_solve_desired_triplets(triplets)

    S_list = list(S)
    S_bar_list = list(S_bar)

    # Extract character matrices corresponding to the partition
    cm = tree.character_matrix
    cm1 = cm.loc[S_list]
    cm2 = cm.loc[S_bar_list]

    # Initialize subtrees
    sub_tree1 = data.CassiopeiaTree(character_matrix=cm1, missing_state_indicator=-1)
    sub_tree2 = data.CassiopeiaTree(character_matrix=cm2, missing_state_indicator=-1)

    # Solve both subtrees using Shared Mutation Joining
    smjsolver = SharedMutationJoiningSolver(
        similarity_function=dissimilarity_functions.hamming_similarity_without_missing
    )

    smjsolver.solve(sub_tree1)
    smjsolver.solve(sub_tree2)

    # Initialize reconstructed tree
    recon_tree = nx.DiGraph()
    connection_node = "root"

    recon_tree.add_node(connection_node)

    # Add subtree 1
    for node in sub_tree1.get_tree_topology().nodes:
        recon_tree.add_node(node)

        for neighbor in sub_tree1.get_tree_topology().neighbors(node):
            recon_tree.add_edge(node, neighbor)

    # Add subtree 2
    for node in sub_tree2.get_tree_topology().nodes:
        recon_tree.add_node(node)

        for neighbor in sub_tree2.get_tree_topology().neighbors(node):
            recon_tree.add_edge(node, neighbor)

    # Connect root to both subtrees
    recon_tree.add_edge(connection_node, list(sub_tree1.get_tree_topology().nodes)[0])
    recon_tree.add_edge(connection_node, list(sub_tree2.get_tree_topology().nodes)[0])

    return recon_tree


def calculate_triplets_correct(
    ground_tree: data.CassiopeiaTree,
    recon: nx.DiGraph
) -> Tuple[int, int]:
    """
    Sample triplets and compute the fraction that are correctly resolved.

    Args:
        ground_tree: Ground truth tree
        recon: Reconstructed tree

    Returns:
        accuracy: Fraction of correctly resolved triplets
    """

    num_correct = 0

    ground = ground_tree.get_tree_topology()

    # Obtain all leaf nodes
    leaves = ground_tree.leaves

    for _ in range(5000):

        # Randomly sample a triplet of leaves
        sampled_triplet = np.random.choice(leaves, 3, replace=False)

        # Determine triplet structure in both trees
        recon_triplet = find_triplet_structure(sampled_triplet, recon)
        ground_triplet = find_triplet_structure(sampled_triplet, ground)

        if recon_triplet == ground_triplet:
            num_correct += 1

    total_triplets = 5000

    accuracy = num_correct / total_triplets if total_triplets > 0 else 0

    return accuracy