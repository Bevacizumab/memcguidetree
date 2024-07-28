# MEMC Guide Tree
-----------------
The minimum evolution with the molecular clock (MEMC) guide tree algorithm.

Guide trees are [hierarchical clustering algorithms](https://en.wikipedia.org/wiki/Hierarchical_clustering) used in [multiple sequence alignments (MSAs)](https://en.wikipedia.org/wiki/Multiple_sequence_alignment).
They direct the order of alignments of input sequences for fast MSAs.

The minimum evolution is related to the minimum total sum of branch lengths.
Therefore, the problem of finding the optimal tree can be mapped into a traveling salesman problem (TSP).

Quantum annealers can provide solutions to TSPs with scalability.
This algorithm can construct guide trees from a given TSP solution with a linear time complexity, O(N).
When constructing guide trees, the molecular clock hypothesis is used.

This guide tree can be used for MAFFT.
Because it requires sequence names in the guide trees to be named as 1, 2, ... n,
we followed that convention in this code.
The guide trees constructed by our algorithm must be converted to the native tree format of MAFFT.
This can be done by the [ruby script](https://mafft.cbrc.jp/alignment/software/newick2mafft.rb) provided by MAFFT.

With simple modifications to accept sequence names, this code is also applicable to various MSA tools.
If you follow the included LICENSE file with a proper citation of the code and our paper, we are happy to let you make any modifications and use this code for other MSA tools.

The tests for this guide tree algorithm can be performed by various MSA benchmarks.
If you want to use [BAliBASE 3.0](https://www.lbgi.fr/balibase/), use [this git repository](https://github.com/robinhundt/bali-score) for calculating scores (Rust compiler is required).



# Example Usage:
--------------
1. memc_guide_tree.py
```python
from memc_guide_tree import MEMCGuideTree

dm_path = 'dm.txt'  # Use actual file path
tsp_path = 'tsp.txt'  # Use actual file path
output_path = 'output.dnd'  # Use actual file path

memc = MEMCGuideTree(dm_path, tsp_path, output_path)
memc.guide_tree_to_file()  # Save the guide tree.
```

memc_guide_tree requires unionfind.py implemented by [Debajyoti Nandi](https://github.com/deehzee/unionfind).
unionfind.py is under MIT license and its copy is provided in this repository as the name file name.


2. tsp_from_dm.py
```python
dm_file_path = 'dm.txt'  # Use actual file path
is_classical = False  # D-Wave Quantum Annealer (LeapHybridSampler).
tsp_order_output_file_path = 'tsp_order.txt'  # Use actual file path

tsp_to_file(dm_file_path, is_classical, tsp_order_output_file_path)
```

From a given distance matrix saved in a file format that supports numpy.loadtxt,
you can find its TSP solution either by D-Wave quantum annealer
```python
is_classical = False
```
or classical algorithm ([python-tsp](https://github.com/fillipe-gsm/python-tsp)).
```python
is_classical = True
```


3. classic_guide_tree.py
```python
algorithm = "upgma"  # Choose a classical algorithm.
dm_path = "dm.txt"  # Use actual file path
output_path = "output.dnd"  # Use actual file path
gt = GuideTree(dm_path, algorithm, output_path)
gt.guide_tree_to_file()  # Save the guide tree.
```

This code provides classical algorithms for guide trees.
All these implementations are based on Biopython with modifications by Youngjun Park(pyoungj09@gmail.com).
- Complete linkage [1]
- Flexible (beta = 0.5) [1]
- Modified UPGMA (in MAFFT) [2]
- Single linkage [1]
- UPGMA [1]
- UPGMC [1]
- Ward [1]
- WPGMA [1]
- WPGMC [1]
are provided in this code.

References to classical algorithms: 
[1] [Milligan, G.W. Ultrametric hierarchical clustering algorithms. Psychometrika 44, 343-346 (1979)](https://doi.org/10.1007/BF02294699)
[2] [Kazutaka Katoh, John Rozewicki, Kazunori D Yamada, MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization, Briefings in Bioinformatics, Volume 20, Issue 4, July 2019, Pages 1160-1166](https://doi.org/10.1093/bib/bbx108)



# Dependencies
--------------
This project requires the following packages:
- [numpy](https://github.com/numpy/numpy)
- [networkx](https://github.com/networkx/networkx)
- [dwave-networkx](https://github.com/dwavesystems/dwave-networkx)
- [dwave-system](https://github.com/dwavesystems/dwave-system)
- [BioPython](https://github.com/biopython/biopython)
- [python-tsp](https://github.com/fillipe-gsm/python-tsp)

I appreciate contributors and maintainers for these valuable packages.



# License
---------
Copyright (c) 2024, Youngjun Park (pyoungj09@gmail.com)

Each code file includes a header explaining the applicable license.
Please see the LICENSE file in the project root for more information.



# Get in Touch
--------------
If you have any questions, don't hesitate to contact 
the author, [Youngjun Park](pyoungj09@gmail.com).



# Citing this work
------------------
(Link will be added).