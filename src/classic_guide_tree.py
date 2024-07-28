# Derived from Biopython's TreeConstruction.py with modifications.
# 
# This source code is licensed under the BSD 3-Clause License. 
# Please see the LICENSE file in the project root for more information.
# 
# Original Biopython Code:
# Copyright (C) 2013 by Yanbo Ye (yeyanbo289@gmail.com)
# 
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Modifications and additional algorithms:
# Copyright (c) 2024, Youngjun Park (pyoungj09@gmail.com)
# All rights reserved.
#
# Biopython UPGMA documentation:
# https://biopython.org/docs/dev/api/Bio.Phylo.TreeConstruction.html
# Biopython UPGMA GitHub:
# https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py#L708



import copy
from Bio.Phylo import BaseTree
import numpy as np
from Bio.Phylo.TreeConstruction import _Matrix
from Bio import Phylo



# Biopython's UPGMA algorithm is in fact WPGMA 
# (Weighted Pair Group Method with Arithmetic Mean), 
# because it simply averages the distances 
# when updating a distance matrix.
# Therefore, to calculate and keep track of the number of nodes
# we make a few modifications to the Biopython's UPGMA algorithm
# with the help of the following three functions written 
# by Youngjun Park (pyoungj09@gmail.com).

def combine_two_elements(lst, index1, index2):
    """
    Combine two specific elements (index1 and index2) 
    in a list (lst) into a sublist.
    Let index1 < index2.
    index2 will be removed from the list.
    
    Return a nested list.
    
    Example Usage:
    --------------
    >>> lst = [[1],[1],[1],[1],[1]]
    >>> lst= combine_two_elements(lst, 0, 1)
    >>> print(lst)
    [[[1], [1]], [1], [1], [1]]
    """
    
    combined = [lst[index1], lst[index2]]
    # Replace the first index with the combined list 
    # and remove the second index (index2)
    lst[index1] = combined
    del lst[index2]
    return lst

def nested_sum(L):
    """
    Return the number of nodes
    Return the sum of all elements in a nested list recursively.
    Because all elements are set to be [1], 
    the sum will be the number of elements in the list.
    
    This function is used in the function count_num_node.
    """
    
    total = 0  # don't use `sum` as a variable name
    for i in L:
        if isinstance(i, list):  # checks if `i` is a list
            total += nested_sum(i)
        else:
            total += i
    return total

def count_num_node(lst, index):
    """
    Count the total number of nodes in each nested list.
    Each nest list represents a subtree.
    When [1] is given, the function returns 1.
    
    Example Usage:
    --------------
    >>> lst = [[[1], [1]], [1], [1], [[[1], [1]], [1]]]
    >>> print(count_num_node(lst, 0))
    2
    >>> print(count_num_node(lst, 1))
    1
    >>> print(count_num_node(lst, 3))
    3
    """

    if len(lst[index]) == 2:
        return (nested_sum(lst[index][0]) 
                + nested_sum(lst[index][1]))
    elif len(lst[index]) == 1:
        return 1



class GuideTree:
    """
    Classical hierarchical clustering algorithms (allowed_algs).
    All these implementations are based on 
    Biopython with modifications.

    Youngjun Park (pyoungj09@gmail.com) added additional 
    classical algorithms except UPGMA, including:
    - Single linkage [1]
    - Complete linkage [1]
    - WPGMA [1]
    - UPGMC [1]
    - WPGMC [1]
    - Ward [1]
    - Flexible [1]
    - Modified UPGMA (in MAFFT) [2]

    All of them were implemented with a time complexity of O(N^3).
    
    References:
    -----------
    [1] Milligan, G.W. Ultrametric hierarchical clustering algorithms. 
        Psychometrika 44, 343-346 (1979). 
        https://doi.org/10.1007/BF02294699
    
    [2] Kazutaka Katoh, John Rozewicki, Kazunori D Yamada, 
        MAFFT online service: multiple sequence alignment, 
        interactive sequence choice and visualization, 
        Briefings in Bioinformatics, 
        Volume 20, Issue 4, July 2019, Pages 1160-1166, 
        https://doi.org/10.1093/bib/bbx108

    Original Biopython code:
    ------------------------
    https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py
        
    Example Usage:
    --------------
    >>> algorithm = "upgma"
    >>> dm_path = "dm.txt"
    >>> output_path = "output.dnd"
    >>> gt = GuideTree(dm_path, algorithm, output_path)
    >>> gt.guide_tree_to_file()  # Save the guide tree.
    """
    
    def __init__(self, 
                 dist_matrix_path: str, 
                 algorithm: str, 
                 output_file_path: str) -> None:
        """
        Initialize the class.

        Args:
            dist_matrix_path (str): a path to the distance matrix.
            algorithm (str): a selected algorithm.
            output_file_path (str): a path to the output file.
        """
        
        self.dist_matrix_path = dist_matrix_path
        self.output_file_path = output_file_path
        
        allowed_algs = ['ward', 'modified_upgma', 
                        'single_linkage', 'complete_linkage', 
                        'upgmc', 'wpgmc', 'upgma', 
                        'wpgma', 'flexible'
        ]

        assert algorithm in allowed_algs, (
            f"Invalid algorithm: {algorithm}. "
            f"Allowed values are: {allowed_algs}"
        )

        self.algorithm = algorithm
        
    def _data_load(self) -> None:
        """
        Load the distance matrix and convert it to a numpy matrix.
        As we will input guide trees to MAFFT,
        the sequence names are 1, 2, 3, ..., n.
        """
        
        self.dist_matrix = np.loadtxt(self.dist_matrix_path)
        
        self.seq_name_list = [
            str(i+1) 
            for i in range(self.dist_matrix.shape[0])
        ]
        
    def _lower_triangular_nested_list(self) -> None:
        """
        Convert the distance matrix (square matrix) 
        to a lower triangular nested list.
        This method is required because _Matrix method in biopython
        requires a lower triangular nested list.
        """

        n = self.dist_matrix.shape[0]
        
        self.lower_tri_dist_list = []
        
        for i in range(n):
            row = []
            for j in range(i + 1):
                row.append(self.dist_matrix[i][j])
            self.lower_tri_dist_list.append(row)
    
    def _tree_construction(self):
        """
        Generate a guide tree of 
        a classical hierarchical clustering algorithm.
        """
        
        self._data_load()
        self._lower_triangular_nested_list()

        # Convert self.lower_tri_dist_list to a _Matrix object.
        # The diagonal elements are all zeros.
        # The following example shows a _Matrix when
        # seq_name_list = ['A', 'B', 'C', 'D', 'E'].
        # A  0.000000
        # B  0.230769    0.000000
        # C  0.384615    0.230769    0.000000
        # D  0.538462    0.538462    0.538462    0.000000
        # E  0.615385    0.384615    0.461538    0.153846    0.000000
        #    A           B           C           D           E

        distance_matrix = _Matrix(self.seq_name_list, 
                                  self.lower_tri_dist_list)

        n = len(self.seq_name_list)

        lst = [[1] for _ in range(n)]
        
        # Make a copy of the distance matrix to be used.
        dm = copy.deepcopy(distance_matrix)
        
        # Initialize terminal clades.
        clades = [BaseTree.Clade(None, name) for name in dm.names]
 
        while len(dm) > 1:
            min_dist = dm[1, 0]
            min_i = 1
            min_j = 0
            # Find the minimum index.
            
            for i in range(1, len(dm)):
                for j in range(0, i):
                    if min_dist >= dm[i, j]:
                        min_dist = dm[i, j]
                        min_i = i
                        min_j = j

            # Create clades.
            clade1 = clades[min_i]
            clade2 = clades[min_j]

            inner_clade = BaseTree.Clade(None, None)
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            
            # Assign branch lengths.
            if clade1.is_terminal():
                clade1.branch_length = min_dist / 2
            else:
                clade1.branch_length = max(
                    0.00001, min_dist / 2 - self._height_of(clade1)
                )
            
            if clade2.is_terminal():
                clade2.branch_length = min_dist / 2
            else:
                clade2.branch_length = max(
                    0.00001, min_dist / 2 - self._height_of(clade2)
                )

            # Update node list.
            clades[min_j] = inner_clade
            
            # Delete min_i.
            del clades[min_i]

            # Rebuild distance matrix.
            # Set the distances of new node at the index of min_j.
            
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    n_min_i = count_num_node(lst, min_i)
                    n_min_j = count_num_node(lst, min_j)
                    n_k = count_num_node(lst, k)
                    
                    if self.algorithm == 'ward':
                        dm[min_j, k] = (
                            (n_min_i + n_k) 
                            / (n_min_i + n_min_j + n_k) 
                            * dm[min_i, k]
                            + (n_min_j + n_k) 
                            / (n_min_i + n_min_j + n_k) 
                            * dm[min_j, k] 
                            - n_k 
                            / (n_min_i + n_min_j + n_k) 
                            * dm[min_i, min_j]
                        )
                    elif self.algorithm == 'modified_upgma':
                        x = 0.1
                        dm[min_j, k] = (
                            x * 0.5 * (dm[min_i, k] + dm[min_j, k]) 
                            + (1-x) * min(dm[min_i, k], dm[min_j, k])
                        )
                    elif self.algorithm == 'single_linkage':
                        dm[min_j, k] = min(
                            dm[min_i, k], dm[min_j, k]
                        )
                    elif self.algorithm == 'complete_linkage':
                        dm[min_j, k] = max(
                            dm[min_i, k], dm[min_j, k]
                        )
                    elif self.algorithm == 'upgmc':
                        dm[min_j, k] = (
                            (n_min_i / (n_min_i + n_min_j)) 
                            * dm[min_i, k] 
                            + (n_min_j / (n_min_i + n_min_j)) 
                            * dm[min_j, k] 
                            - (
                                (n_min_i * n_min_j)
                                / ((n_min_i + n_min_j) ** 2)
                            ) * dm[min_i, min_j]
                        )
                    elif self.algorithm == 'wpgmc':
                        dm[min_j, k] = (
                            0.5 * dm[min_i, k] 
                            + 0.5 * dm[min_j, k] 
                            - 0.25 * dm[min_i, min_j]
                        )
                    elif self.algorithm == 'upgma':
                        dm[min_j, k] = (
                            (n_min_i / (n_min_i + n_min_j))
                            * dm[min_i, k]
                            + (n_min_j / (n_min_i + n_min_j)) 
                            * dm[min_j, k]
                        )
                    elif self.algorithm == 'wpgma':
                        dm[min_j, k] = (
                            0.5 * dm[min_i, k] 
                            + 0.5 * dm[min_j, k]
                        )
                    elif self.algorithm == 'flexible':
                        beta = 0.5
                        dm[min_j, k] = (
                            0.5 * (1 - beta) * dm[min_i, k] 
                            + 0.5 * (1 - beta) * dm[min_j, k] 
                            + beta * dm[min_i, min_j]
                        )
                    else:
                        raise ValueError(
                            f"Invalid algorithm: {self.algorithm}"
                        )

            lst = combine_two_elements(lst, min_i, min_j)
   
            # del dm[min_i] deletes both the row and column 
            # of the index min_i.
            # For example,
            # if m = _Matrix(names=['A', 'B', 'C', 'D'], 
            # matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]]),
            # >>> del m['A']
            # >>> m
            # _Matrix(names=['B', 'C', 'D'], matrix=[[0], [4, 0], [5, 6, 0]])
            del dm[min_i]
            
        inner_clade.branch_length = 0
        return BaseTree.Tree(inner_clade)

    def _height_of(self, clade):
        """
        Calculate the longest height to any terminal clade.
        
        Args:
            clade (BaseTree.Clade): an instance of BaseTree.Clade 
                representing a clade.

        Returns:
            float: The height of the clade.
        """
        
        if clade.is_terminal():
            return clade.branch_length or 0
        else:
            # Recursively calculate the height of the clade.
            return ((clade.branch_length or 0) 
                    + max(self._height_of(c) 
                          for c in clade.clades))

    def guide_tree_to_file(self) -> None:
        """
        Save the guide tree to a file.
        """

        guide_tree = self._tree_construction()
        # output_folder_path ends with a file name (.dnd).
        Phylo.write(guide_tree, self.output_file_path, "newick")