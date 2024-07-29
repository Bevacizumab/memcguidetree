# Copyright (c) 2024, Youngjun Park (pyoungj09@gmail.com)
# 
# This source code is licensed under the BSD 3-Clause License. 
# Please see the LICENSE file in the project root for more information.



from Bio.Phylo import BaseTree
import numpy as np
from Bio import Phylo
from unionfind import UnionFind
import decimal



class MEMCGuideTree:
    """
    This class constructs a minimum evolution 
    with the molecular clock (MEMC) guide tree.
    
    The guide tree is used for MAFFT, so input sequences 
    are named as 1, 2, 3, ..., n.
    However, it can be used for other multiple sequence alignment 
    tools by accepting the input sequences names.

    Example Usage:
    --------------
    >>> from memc_guide_tree import MEMCGuideTree
    >>> dm_path = 'dm.txt'
    >>> tsp_path = 'tsp.txt'
    >>> output_path = 'output.dnd'
    >>> memc = MEMCGuideTree(dm_path, tsp_path, output_path)
    >>> memc.guide_tree_to_file()  # Save the guide tree.
    """
    
    def __init__(self, 
                 dist_matrix_path: str, 
                 tsp_order_path: str,
                 output_file_path: str) -> None:
        """ 
        Initialize the class.

        Args:
            dist_matrix_path (str): a path to the distance matrix file
                (any extension acceptable to numpy.loadtxt).
            tsp_order_path (str): a path to the TSP order file
                (any extension acceptable to numpy.loadtxt).
            output_file_path (str): a path to the output file. 
                The tree is saved as newick format.
        """
        
        self.dm_path = dist_matrix_path
        self.tsp_path = tsp_order_path
        self.output_path = output_file_path

    def _data_load(self) -> None:
        """
        Load the distance matrix and convert it to a numpy matrix.
        Load the tsp order (circular sum order) 
        as a numpy array and convert it to a list.
        """
        
        self.dist_matrix = np.loadtxt(self.dm_path)
        self.tsp_order = np.loadtxt(self.tsp_path, dtype=int)
        
        # TSP order starts from 0 while MAFFT requires 
        # sequnces to be named as 1, 2, 3, ..., n.
        self.tsp_order_seq_name = (
            self.tsp_order + 1
        ).astype(str).tolist()
        self.tsp_order = self.tsp_order.tolist()

    def _counting_sort_for_radix(self, arr, exp):
        """ This method is used in _radix_sort. """
        
        n = len(arr)
        output = [0] * n
        count = [0] * 10

        for i in range(n):
            index = (arr[i][1] // exp) % 10
            count[index] += 1

        for i in range(1, 10):
            count[i] += count[i - 1]

        i = n - 1
        while i >= 0:
            index = (arr[i][1] // exp) % 10
            output[count[index] - 1] = arr[i]
            count[index] -= 1
            i -= 1

        for i in range(n):
            arr[i] = output[i]

    def _radix_sort(self, arr):
        """ 
        Radix sort executed by _sorting_dm.

        Args:
            arr (list): all distances must be scaled to integers.
        """
        
        max1 = max(arr, key=lambda x: x[1])[1]
        
        # Starting from the least significant digit, 
        # sort the distances.
        exp = 1
        
        while max1 / exp > 1:
            self._counting_sort_for_radix(arr, exp)
            
            # Move to the next significant digit.
            exp *= 10

    def _sorting_dm(self, dm, order):
        """ 
        Radix sort the distances between the TSP order.
            This method uses _counting_sort_for_radix and _radix_sort.
            
        Args:
            dm (numpy.ndarray): Loaded distance matrix.
            order (list): Loaded TSP order.
        """
        
        # Make a list of tuples (index, distance).
        dm_lst = [
            (i, dm[order[i]][order[(i + 1) % self._n]]) 
            for i in range(self._n)
        ]
        
        # Scale the distances to integers.
        scale_factor = 10 ** len(str(dm[0][1]).split('.')[1])
        scaled_arr = [
            (
                element[0], 
                int(decimal.Decimal(str(element[1])) * scale_factor)
            )
            for element in dm_lst
        ]
        
        # Execute the radix sort.
        self._radix_sort(scaled_arr)
        
        # Rescale the distances to floats.
        self.sorted_dm = [
            (element[0], element[1] / scale_factor) 
            for element in scaled_arr
        ]

    def _memc_guidetree(self):
        """
        Construct the MEMC guide tree.

        This method (_memc_guidetree) is a new algorithm and its code 
        implementation was inspired by Biopython's UPGMA class.
        Specifically, it uses BaseTree.Clade in a similar way:
        https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py#L708
        """
    
        # Step 1. Load the distance matrix and the TSP order.
        # self.dist_matrix, self.tsp_order, self.tsp_order_seq_name.
        self._data_load()
        
        # Number of sequences.
        self._n = len(self.tsp_order_seq_name)
        
        # Step 2. Sort the distances between the TSP order 
        # by radix sort.
        self._sorting_dm(self.dist_matrix, self.tsp_order)

        # Step 3. Initialize the UnionFind data structure, 
        # clades, and height.
        uf = UnionFind(self.tsp_order_seq_name)
        clades = [
            BaseTree.Clade(None, self.tsp_order_seq_name[i])
            for i in range(self._n)
        ]
        height = np.zeros(self._n, dtype=float)
        
        # Step 4. Construct the guide tree.
        # Repeat until there is only one clade left.
        while len(self.sorted_dm) > 1:
            
            # (index)-th and (index+1)-th sequences are united.
            # The distance between them is dist.
            index, dist = self.sorted_dm.pop(0)
            
            # Because the instance of UnionFind was initialized 
            # by self.tsp_order_seq_name,
            # convert the index to the name of the sequence.
            tsp_1 = self.tsp_order_seq_name[index]
            tsp_2 = self.tsp_order_seq_name[(index + 1) % self._n]
            
            # We use weighted-quick-union-with-path-compression.
            # Git: https://github.com/deehzee/unionfind
            # Ref: http://algs4.cs.princeton.edu/lectures/
            uf_root_1 = uf.find(tsp_1)
            uf_root_2 = uf.find(tsp_2)
            
            uf.union(tsp_1, tsp_2)
            united_root = uf.find(tsp_1)
            
            # deleted_root is the root that was merged.
            if united_root == uf_root_1:
                deleted_root = uf_root_2
            else:
                deleted_root = uf_root_1

            clade1 = clades[uf_root_1]
            clade2 = clades[uf_root_2]
            
            # Make a new inner clade.
            inner_clade = BaseTree.Clade(None, None)
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            
            # Assign branch length.
            if clade1.is_terminal():
                clade1.branch_length = dist / 2
            else:
                # Ensure non-zero branch length.
                clade1.branch_length = max(
                    0.00001, dist / 2 - height[uf_root_1]
                )
            
            if clade2.is_terminal():
                clade2.branch_length = dist / 2
            else:
                # Ensure non-zero branch length.
                clade2.branch_length = max(
                    0.00001, dist / 2 - height[uf_root_2]
                )
            
            clades[united_root] = inner_clade
            clades[deleted_root] = None
            
            if height[deleted_root] > 0:
                height[deleted_root] = 0.0

            # Save the height of the united clade.
            height[united_root] = dist / 2

        # Set the branch length of the root clade to 0.
        inner_clade.branch_length = 0
        
        # Return the guide tree.
        tree = BaseTree.Tree(inner_clade)
        return tree
    
    def guide_tree_to_file(self) -> None:
        """
        Generate the memc guide tree and write it to a file.
        """
        tree = self._memc_guidetree()
        
        # output_folder_path ends with a file name 
        # (for example, .dnd).
        Phylo.write(tree, self.output_path, "newick")