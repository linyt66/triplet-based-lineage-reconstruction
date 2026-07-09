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

def percolation_solve(tree: data.CassiopeiaTree) -> data.CassiopeiaTree:
    """
    Reconstructs a tree using the Percolation (Threshold) Algorithm.

    Parameters
    ----------
    tree : CassiopeiaTree
        Input tree containing only a character matrix.

    Returns
    -------
    CassiopeiaTree
        Reconstructed tree topology.
    """
    cm = tree.character_matrix.astype(int)

    recon_tree = data.CassiopeiaTree(
        character_matrix=cm,
        missing_state_indicator=-1
    )

    solver = PercolationSolver(
        joining_solver=VanillaGreedySolver()
    )
    solver.solve(recon_tree)

    return recon_tree

def shared_mutation_solve(tree: data.CassiopeiaTree) -> data.CassiopeiaTree:
    """
    Reconstructs a tree using Shared Mutation Joining.

    Parameters
    ----------
    tree : CassiopeiaTree
        Input tree containing only a character matrix.

    Returns
    -------
    CassiopeiaTree
        Reconstructed tree topology.
    """
    cm = tree.character_matrix.astype(int)

    recon_tree = data.CassiopeiaTree(
        character_matrix=cm,
        missing_state_indicator=-1
    )

    solver = SharedMutationJoiningSolver(
        similarity_function=dissimilarity_functions.hamming_similarity_without_missing
    )
    solver.solve(recon_tree)

    return recon_tree

def robinson_foulds_zero(
    ground_tree: data.CassiopeiaTree,
    recon_tree: data.CassiopeiaTree
) -> bool:
    """
    Checks whether two trees have zero Robinson–Foulds distance.

    Returns
    -------
    bool
        True if RF distance equals zero.
    """
    rf, _ = critique.compare.robinson_foulds(
        ground_tree, recon_tree
    )
    return rf == 0

def find_triplet_structure(
    triplet: tuple[str, str, str],
    tree: nx.DiGraph
) -> str:
    """
    Determines the ingroup of a triplet in a rooted tree.

    Parameters
    ----------
    triplet : tuple
        Three leaf labels (a, b, c).
    tree : nx.DiGraph
        Tree topology.

    Returns
    -------
    str
        'ab', 'ac', 'bc', or '-' if unresolved.
    """
    a, b, c = triplet

    ancestors = {
        a: set(nx.ancestors(tree, a)),
        b: set(nx.ancestors(tree, b)),
        c: set(nx.ancestors(tree, c)),
    }

    shared = {
        "ab": len(ancestors[a] & ancestors[b]),
        "ac": len(ancestors[a] & ancestors[c]),
        "bc": len(ancestors[b] & ancestors[c]),
    }

    best = max(shared, key=shared.get)
    values = list(shared.values())

    return best if values.count(shared[best]) == 1 else "-"

def triplet_accuracy(
    ground_tree: data.CassiopeiaTree,
    recon_tree: nx.DiGraph,
    n_samples: int = 5000
) -> float:
    """
    Estimates triplet accuracy by random sampling.

    Parameters
    ----------
    ground_tree : CassiopeiaTree
        Ground truth tree.
    recon_tree : nx.DiGraph
        Reconstructed tree.
    n_samples : int
        Number of triplets sampled.

    Returns
    -------
    float
        Proportion of correctly resolved triplets.
    """
    ground = ground_tree.get_tree_topology()
    leaves = list(ground_tree.leaves)

    correct = 0
    for _ in range(n_samples):
        triplet = np.random.choice(leaves, 3, replace=False)
        correct += (
            find_triplet_structure(triplet, recon_tree)
            == find_triplet_structure(triplet, ground)
        )

    return correct / n_samples

def hamming_distance(x, y) -> int:
    """Computes Hamming distance between two vectors."""
    return np.sum(x != y)

def infer_triplets_from_mutations(
    tree: data.CassiopeiaTree
) -> list[tuple[str, str, str]]:
    """
    Infers triplets based on minimum pairwise Hamming distance.

    Returns
    -------
    list
        Triplets of the form (ingroup1, ingroup2, outgroup).
    """
    cm = tree.character_matrix
    triplets = []

    for a, b, c in itertools.combinations(tree.leaves, 3):
        dists = {
            (a, b): hamming_distance(cm.loc[a], cm.loc[b]),
            (a, c): hamming_distance(cm.loc[a], cm.loc[c]),
            (b, c): hamming_distance(cm.loc[b], cm.loc[c]),
        }

        (x, y), _ = min(dists.items(), key=lambda t: t[1])
        outgroup = ({a, b, c} - {x, y}).pop()

        triplets.append((x, y, outgroup))

    return triplets

def construct_triplet_graph(
    triplets: list[tuple[str, str, str]]
) -> nx.Graph:
    """
    Constructs a weighted graph encoding triplet constraints.
    """
    G = nx.Graph()

    for a, b, c in triplets:
        G.add_edge(a, b, weight=G.get_edge_data(a, b, {"weight": 0})["weight"] - 2)
        G.add_edge(a, c, weight=G.get_edge_data(a, c, {"weight": 0})["weight"] + 1)
        G.add_edge(b, c, weight=G.get_edge_data(b, c, {"weight": 0})["weight"] + 1)

    return G

def gw_partition(triplets):
    """
    Applies Goemans–Williamson SDP to partition leaves.
    """
    G = construct_triplet_graph(triplets)
    cut = cvxgr.algorithms.goemans_williamson_weighted(G)
    return set(cut.left), set(cut.right)

def build_tree_from_triplets(
    tree: data.CassiopeiaTree,
    triplets
) -> nx.DiGraph:
    """
    Builds a hierarchical tree using GW partition + recursive SMJ.
    """
    S, S_bar = gw_partition(triplets)
    cm = tree.character_matrix

    subtrees = []
    for subset in (S, S_bar):
        sub_cm = cm.loc[list(subset)]
        sub_tree = data.CassiopeiaTree(
            character_matrix=sub_cm,
            missing_state_indicator=-1
        )
        shared_mutation_solve(sub_tree)
        subtrees.append(sub_tree.get_tree_topology())

    recon = nx.DiGraph()
    root = "root"
    recon.add_node(root)

    for subtree in subtrees:
        recon.add_nodes_from(subtree.nodes)
        recon.add_edges_from(subtree.edges)
        recon.add_edge(root, list(subtree.nodes)[0])

    return recon
