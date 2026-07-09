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

def percolation_solve(tree: data.CassiopeiaTree) -> data.CassiopeiaTree:
    """Solves a tree using the Threshold Algorithm from a character matrix.

    Args:
        tree: A CassiopeiaTree object without the `tree` instance populated

    Returns:
        A CassiopeiaTree object containing the solved tree topology
    """

    # Extract the character matrix from the CassiopeiaTree container object
    cm = tree.character_matrix
    cm = cm.astype(int)

    # Instantiate a new CassiopeiaTree with no tree initialized
    recon_tree = data.CassiopeiaTree(character_matrix = cm, missing_state_indicator = -1)

    # Instantiate a Threshold Algorithm solver object and solve
    psolver = PercolationSolver(joining_solver = VanillaGreedySolver())
    psolver.solve(recon_tree)

    return recon_tree

def shared_mutation_solve(tree: data.CassiopeiaTree) -> data.CassiopeiaTree:
    """Solves a tree using the Threshold Algorithm from a character matrix.

    Args:
        tree: A CassiopeiaTree object without the `tree` instance populated

    Returns:
        A CassiopeiaTree object containing the solved tree topology
    """
    # Extract the character matrix from the CassiopeiaTree container object
    cm = tree.character_matrix
    cm = cm.astype(int)

    # Instantiate a new CassiopeiaTree with no tree initialized
    recon_tree = data.CassiopeiaTree(character_matrix = cm, missing_state_indicator = -1)

    # Instantiate a Bottom-Up Algorithm solver object and solve
    smjsolver = SharedMutationJoiningSolver(similarity_function = dissimilarity_functions.hamming_similarity_without_missing)
    smjsolver.solve(recon_tree)

    return recon_tree
#求解树的算法，使用了共享突变位点求解
def robinson_foulds_score(ground_tree: data.CassiopeiaTree, recon_tree: data.CassiopeiaTree) -> bool:
    """Determines if the Robinson-Foulds Distance of two tree topologies is 0.

    Args:
        ground_tree: The ground truth tree topology
        recon_tree: The reconstructed tree topology

    Returns:
        If the two topologies have a Robinson-Foulds Distance of 0
    """
    rf, rf_max = critique.compare.robinson_foulds(ground_tree, recon_tree)

    return rf == 0

def find_triplet_structure(triplet: Tuple[str, str, str], tree: nx.DiGraph) -> str:
    """Finds the ingroup of a triplet on a networkx tree topology.

    Args:
        triplet: A tuple containing the three leaf nodes forming the triplet
        tree: The tree topology, stored as a networkx object

    Returns:
        A string determining the ingroup of the triplet
    """
    # Extract the nodes from the tuple
    a, b, c = triplet[0], triplet[1], triplet[2]

    # Find the list of ancestors on the tree for each node
    a_ancestors = [node for node in nx.ancestors(tree, a)]
    b_ancestors = [node for node in nx.ancestors(tree, b)]
    c_ancestors = [node for node in nx.ancestors(tree, c)]

    # Calculate the number of ancestors shared between each leaf node
    ab_common = len(set(a_ancestors) & set(b_ancestors))
    ac_common = len(set(a_ancestors) & set(c_ancestors))
    bc_common = len(set(b_ancestors) & set(c_ancestors))

    # The nodes that share the most number of ancestors are the ingroup,
    # if it is a tie, there is no ingroup
    structure = "-"
    if ab_common > bc_common and ab_common > ac_common:
        structure = "ab"
    elif ac_common > bc_common and ac_common > ab_common:
        structure = "ac"
    elif bc_common > ab_common and bc_common > ac_common:
        structure = "bc"

    return structure

def calculate_triplets_correct(ground_tree: data.CassiopeiaTree, recon_tree: data.CassiopeiaTree, sample: int = 1000) -> Tuple[int, int]:
    """Samples triplets and returns the number that are correctly resolved.

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
    recon = recon_tree.get_tree_topology()
    ground = ground_tree.get_tree_topology()

    leaves = ground.leaves()

    for _ in range(sample):
        # Sample a triplet from the leaves of the ground truth
        sampled_triplet = np.random.choice(leaves, 3, replace = False)
        # Find the ingroup for both tree topologies
        recon_triplet = find_triplet_structure(sampled_triplet, recon)
        ground_triplet = find_triplet_structure(sampled_triplet, ground)
        # If the ingroup is the same, add one to the number of correct triplets
        num_correct += int(recon_triplet == ground_triplet)

    return num_correct, sample

def triplets_score(ground_tree: data.CassiopeiaTree, recon_tree: data.CassiopeiaTree) -> bool:
    """Determines if the proportion of triplets correctly resolved is >= 95%.

    Args:
        ground_tree: The ground truth tree topology
        recon_tree: The reconstructed tree topology

    Returns:
        If the topologies have the same ingroup on >= 95% of sampled triplets
    """
    correct_triplets, total_triplets = calculate_triplets_correct(ground_tree, recon_tree)

    return correct_triplets/total_triplets >= 0.95
    
def depth_isomorphism(d: float, ground_tree: data.CassiopeiaTree, recon_tree: data.CassiopeiaTree) -> bool:
    """Determines if two trees have the same splits up to depth d.

    Determining if two trees have the same splits up to depth d (are depth 
    isomorphic at d) is the same as saying that all triplets whose LCAs
    occur before depth d will be resolved correctly. The reconstruction 
    criterion is the same as setting \ell^* equal to the minimum edge length
    and d^* to d, as are done in the simulations varying d^*.

    Args:
        d: The depth to which to check if the splits are the same for two trees
        ground_tree: The ground truth tree topology
        recon_tree: The reconstructed tree topology

    Returns:
        If the topologies have the same splits up to depth d
    """

    tree_g = copy.copy(ground_tree)
    tree_r = copy.copy(recon_tree)

    tree_g.collapse_unifurcations()
    tree_r.collapse_unifurcations()
    
    def check_split(node_g, tree_g, node_r, tree_r, d, flag):
        groups_g = {}
        for i in tree_g.children(node_g):
            groups_g[i] = set(tree_g.leaves_in_subtree(i))
        groups_r = {}
        for i in tree_r.children(node_r):
            groups_r[i] = set(tree_r.leaves_in_subtree(i))
        node_pairs = []
        for i in groups_g:
            for j in groups_r:
                if groups_g[i] == groups_r[j]:
                    node_pairs.append((i, j))
        if len(node_pairs) == 0:
            flag[0] = False
        for i in node_pairs:
            if tree_g.get_attribute(i[0], "time") <= d and flag[0] and len(tree_g.leaves_in_subtree(i[0])) > 1:
                check_split(i[0], tree_g, i[1], tree_r, d, flag)
    
    correct_flag = [True]            
    check_split(tree_g.root, tree_g, tree_r.root, tree_r, d, correct_flag)

    return correct_flag[0]

def complete_binary_topology_sim(depth: int) -> data.CassiopeiaTree:
    """Generates a complete binary topology.

    Used in the validation of k experiments where multiple mutation datasets
    are generated for each topology.
    
    Args:
        depth: The depth of the tree (including the unary implicit root edge)

    Returns:
        A complete binary topology
    """

    # Creating the binary networkx tree   
    topology = nx.balanced_tree(2, (depth - 1), create_using = nx.DiGraph)
    topology = nx.relabel_nodes(topology, dict(zip(list(topology.nodes), [i + 1 for i in topology.nodes])))

    # Appending the implicit root edge
    topology.add_edge(0, 1)

    # Storing the networkx tree in a CassiopeiaTree container object
    complete_binary = data.CassiopeiaTree(tree = topology)

    # Relabling the leaves
    leaves = complete_binary.leaves
    complete_binary.relabel_nodes(dict(zip(leaves, ['c' + str(i) for i in leaves])))

    # Normalizing the times and edge lengths by the total time of the tree
    time_dict = {}
    for i in complete_binary.nodes:
        time_dict[i] = complete_binary.get_time(i)/depth
    complete_binary.set_times(time_dict)

    return complete_binary

def exponential_plus_c_topology_sim(depth: int) -> data.CassiopeiaTree:
    """Generates an asynchronous topology.

    Used in the validation of k experiments where multiple mutation datasets
    are generated for each topology.

    Args:
        depth: The depth of the tree (including the unary implicit root edge)

    Returns:
        An asynchronous topology
    
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
    
    return topology

def overlay_mut_data(topology: data.CassiopeiaTree, k_cand: int, q_dist: Dict[int, float], lamb: float) -> data.CassiopeiaTree:
    """Overlays lineage tracing mutations over a provided topology.

    Used in the validation of k experiments where multiple mutation datasets
    are generated for each topology.

    Args:
        topology: The CassiopeiaTree object containing the tree topology
        k_cand: The number of characters
        q_dist: The state distribution
        lamb: The value of lambda

    Returns:
        A CassiopeiaTree object containing the topology as well as the
        lineage tracing mutations
    """

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

def find_recon_triplets(tree:data.CassiopeiaTree):
#生成以共享突变位点为依据的重建树
   cm = tree.character_matrix  # 获取突变信息矩阵
   results = []

    # 遍历所有三元组组合
   for a, b, c in combinations(tree.leaves, 3):
        # 提取对应突变信息
        mutations_a = cm.loc[a]
        mutations_b = cm.loc[b]
        mutations_c = cm.loc[c]

        # 计算汉明距离
        distance_ab = hamming_distance(mutations_a, mutations_b)
        distance_ac = hamming_distance(mutations_a, mutations_c)
        distance_bc = hamming_distance(mutations_b, mutations_c)

        # 找到汉明距离最小的两个元素
        min_distance = min(distance_ab, distance_ac, distance_bc)

        if min_distance == distance_ab:
            # a和b是最小距离的那对，c就是外群
            outer_group = c
            results.append((a, b, outer_group))
        elif min_distance == distance_ac:
            # a和c是最小距离的那对，b就是外群
            outer_group = b
            results.append((a,c, outer_group))
        else:
            # b和c是最小距离的那对，a就是外群
            outer_group = a
            results.append((b, c, outer_group))

        # 记录结果
   return results

def find_triplet_structure(triplet: Tuple[str, str, str], tree: nx.DiGraph) -> str:
    """Finds the ingroup of a triplet on a networkx tree topology.

    Args:
        triplet: A tuple containing the three leaf nodes forming the triplet
        tree: The tree topology, stored as a networkx object

    Returns:
        A string determining the ingroup of the triplet
    """
    # Extract the nodes from the tuple
    a, b, c = triplet[0], triplet[1], triplet[2]

    # Find the list of ancestors on the tree for each node
    a_ancestors = [node for node in nx.ancestors(tree, a)]
    b_ancestors = [node for node in nx.ancestors(tree, b)]
    c_ancestors = [node for node in nx.ancestors(tree, c)]

    # Calculate the number of ancestors shared between each leaf node
    ab_common = len(set(a_ancestors) & set(b_ancestors))
    ac_common = len(set(a_ancestors) & set(c_ancestors))
    bc_common = len(set(b_ancestors) & set(c_ancestors))

    # The nodes that share the most number of ancestors are the ingroup,
    # if it is a tie, there is no ingroup
    structure = "-"
    if ab_common > bc_common and ab_common > ac_common:
        structure = "ab"
    elif ac_common > bc_common and ac_common > ab_common:
        structure = "ac"
    elif bc_common > ab_common and bc_common > ac_common:
        structure = "bc"
    return structure

def hamming_distance(vec1, vec2):
    """计算两个行向量之间的汉明距离"""
    return np.sum(vec1 != vec2)

def find_recon_triplets(tree:data.CassiopeiaTree):
#生成以共享突变位点为依据的重建树
   cm = tree.character_matrix  # 获取突变信息矩阵
   results = []

    # 遍历所有三元组组合
   for a, b, c in combinations(tree.leaves, 3):
        # 提取对应突变信息
        mutations_a = cm.loc[a]
        mutations_b = cm.loc[b]
        mutations_c = cm.loc[c]

        # 计算汉明距离
        distance_ab = hamming_distance(mutations_a, mutations_b)
        distance_ac = hamming_distance(mutations_a, mutations_c)
        distance_bc = hamming_distance(mutations_b, mutations_c)

        # 找到汉明距离最小的两个元素
        min_distance = min(distance_ab, distance_ac, distance_bc)

        if min_distance == distance_ab:
            # a和b是最小距离的那对，c就是外群
            outer_group = c
            results.append((a, b, outer_group))
        elif min_distance == distance_ac:
            # a和c是最小距离的那对，b就是外群
            outer_group = b
            results.append((a,c, outer_group))
        else:
            # b和c是最小距离的那对，a就是外群
            outer_group = a
            results.append((b, c, outer_group))

        # 记录结果
   return results

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