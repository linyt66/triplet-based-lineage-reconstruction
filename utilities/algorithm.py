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
    Construct a weighted graph based on the triplets.

    Args:
    - triplets: List of triplets, where each triplet is a tuple of three node values.

    Returns:
    - G: Weighted graph constructed from the triplets using NetworkX.
    """
    G = nx.MultiGraph()
    for triplet in triplets:
        a, b, c = triplet
        # Add edges with weights according to the constraints for TD triplets
        G.add_edge(a, b, weight=-2)
        G.add_edge(a, c, weight=1)
        G.add_edge(b, c, weight=1)

    # 转换为简单图并累加权重
    G_simple = nx.Graph()
    for u, v, data in G.edges(data=True):
        if G_simple.has_edge(u, v):
            G_simple[u][v]['weight'] += data['weight']  # 累加权重
        else:
            G_simple.add_edge(u, v, weight=data['weight'])  # 添加边及其初始权重

    return G_simple

def goemans_williamson_solve_desired_triplets(triplets):
    G = construct_graph(triplets)
    sdp_cut = cvxgr.algorithms.goemans_williamson_weighted(G)
    best_partition=[set(), set()]
    best_partition[0]=sdp_cut.left
    best_partition[1]=sdp_cut.right
    S=set(best_partition[0])
    S_bar=set(best_partition[1])
    obeyed_triplets = []
    disobeyed_triplets = []
    for triplet in triplets:
        a, b, c = triplet
        if (a in S and b in S and c not in S) or (a not in S and b not in S and c in S):
            obeyed_triplets.append(triplet)
        else:
            disobeyed_triplets.append(triplet)
    cut_value=sdp_cut.evaluate_cut_size(G)
  #  accuarcy1=len(obeyed_triplets) / len(errors_triplets)
   # accuarcy=(1/3*len(errors_triplets)+1/3*cut_value)/len(errors_triplets)
    return S, S_bar


def build_tree_from_triplet_partition(tree: data.CassiopeiaTree, triplets) -> nx.DiGraph:
    S, S_bar = goemans_williamson_solve_desired_triplets(triplets)
    S_list = list(S)
    S_bar_list = list(S_bar)

    # Extract the character matrices from the partition
    cm = tree.character_matrix
    cm1 = cm.loc[S_list]  # Use list to index
    cm2 = cm.loc[S_bar_list]
    
    # Initialize subtrees
    sub_tree1 = data.CassiopeiaTree(character_matrix=cm1, missing_state_indicator=-1)
    sub_tree2 = data.CassiopeiaTree(character_matrix=cm2, missing_state_indicator=-1)
    
    # Solve for both subtrees
    smjsolver = SharedMutationJoiningSolver(similarity_function=dissimilarity_functions.hamming_similarity_without_missing)
    smjsolver.solve(sub_tree1)
    smjsolver.solve(sub_tree2)

    # Create a directed graph for the reconstruction tree
    recon_tree = nx.DiGraph()
    connection_node = "root"  # Define the connection node name

    # Add the connection node
    recon_tree.add_node(connection_node)

    # Add sub_tree1 and sub_tree2 to the reconstruction tree
    for node in sub_tree1.get_tree_topology().nodes:
        recon_tree.add_node(node)
        for neighbor in sub_tree1.get_tree_topology().neighbors(node):
            recon_tree.add_edge(node, neighbor)

    for node in sub_tree2.get_tree_topology().nodes:
        recon_tree.add_node(node)
        for neighbor in sub_tree2.get_tree_topology().neighbors(node):
            recon_tree.add_edge(node, neighbor)

    # Connect the root to both subtrees
    recon_tree.add_edge(connection_node, list(sub_tree1.get_tree_topology().nodes)[0])  # Connect to first node of sub_tree1
    recon_tree.add_edge(connection_node, list(sub_tree2.get_tree_topology().nodes)[0])  # Connect to first node of sub_tree2

    return recon_tree  # Return the reconstruction tree as a DiGraph

def calculate_triplets_correct(ground_tree: data.CassiopeiaTree, recon: nx.DiGraph) -> Tuple[int, int]:
    """Samples triplets and returns the number that are corrp'pppectly resolved.

    Args:
        ground_tree: The ground truth tree topology
        recon_tree: The reconstructed tree topology
        sample: The number of triplets to sample

    Returns:
        A tuple containing the number of correct triplets and the number sampled
    """

    # Track the number of correct triplets
    num_correct = 0

    # Extract tree topology as networkx objects
   # recon = recon_tree.get_tree_topology()
    ground = ground_tree.get_tree_topology()
    leaves = ground_tree.leaves  #查找所有叶子
    all_triplets = itertools.combinations(leaves, 3)
    
    for _ in range(5000):
        # Sample a triplet from the leaves of the ground truth
        sampled_triplet = np.random.choice(leaves, 3, replace = False)
        # Find the ingroup for both tree topologies
        recon_triplet = find_triplet_structure(sampled_triplet, recon)
        ground_triplet = find_triplet_structure(sampled_triplet, ground)
        # If the ingroup is the same, add one to the number of correct triplets
        num_correct += int(recon_triplet == ground_triplet)
    
    total_triplets = 5000
    '''
    for sampled_triplet in all_triplets:
        total_triplets += 1

        # Find the ingroup for both tree topologies
        recon_triplet = find_triplet_structure(sampled_triplet, recon)
        ground_triplet = find_triplet_structure(sampled_triplet, ground)

        # If the ingroup is the same, increment the count of correct triplets
        if recon_triplet == ground_triplet:
            num_correct += 1
    '''
    # Calculate accuracy
    accuracy = num_correct / total_triplets if total_triplets > 0 else 0
    
    return  accuracy

