# Copyright (c) 2024, Youngjun Park (pyoungj09@gmail.com)
#
# This source code is licensed under the BSD 3-Clause License.
# Please see the LICENSE file in the project root for more information.

"""
TSP solution (minimum evolution circular order, MECO)
from distance matrices.
"""

import networkx as nx
import dwave_networkx.algorithms.tsp as dnx
import numpy as np
from dwave.system import LeapHybridSampler
# classical TSP solver
from python_tsp.exact import solve_tsp_dynamic_programming


def create_edge_weight_list(dist_mat: np.ndarray) -> list:
    """
    Create edge weight list from distance matrix.

    Args:
        dist_mat (np.ndarray): distance matrix.

    Returns:
        list: a list of tuples (i, j, dist_mat[i, j]).
    """

    n = len(dist_mat)
    edge_weight_list = []

    for i in range(n):
        for j in range(i + 1, n):
            edge_weight_list.append((i, j, dist_mat[i, j]))
    return edge_weight_list


def tsp(dm_file_path: str, is_classical: bool) -> list:
    """
    TSP algorithm.

    Example Distance Matrix:
    ------------------------
    0.00000e+00 1.63044e+00 1.56236e+00 1.71959e+00 \\
    1.63044e+00 0.00000e+00 1.56252e+00 1.40315e+00 \\
    1.56236e+00 1.56252e+00 0.00000e+00 1.64347e+00 \\
    1.71959e+00 1.40315e+00 1.64347e+00 0.00000e+00

    Args:
        dm_file_path (str): a path to the distance matrix file.
        is_classical (bool): True for classical, False for quantum.

    Returns:
        list: a TSP order (Minimum evolution circular order).
    """

    # Load distance matrix.
    dm = np.loadtxt(dm_file_path)

    # Create a list of tuples (i, j, dist_mat[i, j]).
    edge_wght_lst = create_edge_weight_list(dm)

    # Initialize a networkx graph and add edge_wght_lst.
    G = nx.Graph()
    G.add_weighted_edges_from(edge_wght_lst)

    if is_classical:
        # Classical TSP algorithm.
        tsp_order, _ = solve_tsp_dynamic_programming(dm)

        return tsp_order
    else:
        # Quantum annealing (D-Wave) TSP algorithm.
        # We recommend using the LeapHybridSampler
        # for large problems. The capacity is determined by
        # a solver's maximum number of binary variables.
        tsp_order = dnx.traveling_salesperson(
            G, LeapHybridSampler(), start=0
        )

        return tsp_order


def tsp_to_file(dm_file_path: str,
                is_classical: bool,
                tsp_order_output_file_path: str):
    """
    Execute the TSP algorithm and save the TSP order to a file.

    Args:
        dm_file_path (str): a path to the distance matrix file.
        is_classical (bool): a flag indicating whether to use
            a classical or a quantum TSP solver.
        tsp_order_output_file_path (str): a path to the output file
            where the TSP order will be saved.

    Example Usage:
    --------------
    >>> from tsp_from_dm import *
    >>> dm_file_path = 'dm.txt'
    >>> is_classical = False  # D-Wave Quantum Annealer.
    >>> tsp_order_output_file_path = 'tsp_order.txt'
    >>> tsp_to_file(dm_file_path, is_classical,
    ...             tsp_order_output_file_path)
    """

    # Execute the TSP algorithm.
    tsp_order = tsp(dm_file_path, is_classical)

    # Save the TSP order to a file.
    np.savetxt(tsp_order_output_file_path, tsp_order, fmt='%d')
