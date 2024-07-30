# MEMC Guide Tree
The Minimum Evolution with the Molecular Clock (MEMC) guide tree algorithm.

Guide trees are [hierarchical clustering algorithms](https://en.wikipedia.org/wiki/Hierarchical_clustering) used in [multiple sequence alignments (MSAs)](https://en.wikipedia.org/wiki/Multiple_sequence_alignment).
They direct the order of alignments of input sequences in MSAs.

The minimum evolution is related to the minimum total sum of branch lengths. The problem of finding the optimal tree corresponding to the minimum evolution can be mapped into a traveling salesman problem (TSP).

Quantum annealers can provide solutions to TSPs with scalability. Our algorithm, implemented in `src/memc_guide_tree.py`, can construct a guide tree from a given TSP solution with a linear time complexity, O(N).
When constructing guide trees, the molecular clock hypothesis is used.

We designed our algorithm to produce guide trees that can be used for [MAFFT](https://mafft.cbrc.jp/alignment/software/).
Because it requires sequence names in the guide trees to be named as 1, 2, ... n, we followed that convention in this code.
The guide trees constructed by our algorithm must be converted to the native tree format of MAFFT.
This can be done using the [ruby script](https://mafft.cbrc.jp/alignment/software/newick2mafft.rb) provided by MAFFT.

With simple modifications to accept sequence names, this code is also applicable to various MSA tools.
If you follow the included LICENSE file with proper citation of the code and our paper, we are happy to let you make any necessary modifications and use this code for other MSA tools.

The tests for this guide tree algorithm can be performed by various MSA benchmarks.
If you want to use [BAliBASE 3.0](https://www.lbgi.fr/balibase/), use [this git repository](https://github.com/robinhundt/bali-score) for calculating scores (Rust compiler is required).



# Example Usage:

This project constructs guide trees for [MAFFT](https://mafft.cbrc.jp/alignment/software/).



## Step 1: Get a Distance Matrix from an MSA Program

To generate a distance matrix using MAFFT, use the following commands:

For FFT-NS-i:
```bash
mafft --retree 0 --maxiterate 0 --distout input_file
```
For L-INS-i:
```bash
mafft --localpair --maxiterate 0 --distout input_file
```
Note that there should be no output file directory in the command.
This command generates a distance matrix in a `.hat2` file in the same directory as the input file.

A typical output format:
```
    1
    (# of seqs)
 (some float)
   1. =(seq1_name)
   2. =(seq2_name)
   3. =(seq3_name)
   4. =(seq4_name)
 d(1,2) d(1,3) d(1,4)
 d(2,3) d(2,4)
 d(3,4)
```
where `d(i,j)` represents the distance between the i-th and the j-th sequences.
The default option outputs distances to the three decimal points.
Modify the `WriteFloatHat2_pointer_halfmtx` function in `io.c` to print more decimal points, at least five or six, for more accurate guide trees.

`memc_guide_tree.py` and `classic_guide_tree.py` expect a distance matrix saved as, for example:
```
0.00000 d(1,2)  d(1,3)  d(1,4)
d(1,2)  0.00000 d(2,3)  d(2,4)
d(1,3)  d(2,3)  0.00000 d(3,4)
d(1,4)  d(2,4)  d(3,4)  0.00000
```



## Step 2. Solve a TSP

Use `tsp_from_dm.py` to solve the TSP:
```python
from tsp_from_dm.py import *

dm_file_path = 'dm.txt'  # Use actual file path
is_classical = False  # D-Wave Quantum Annealer (LeapHybridSampler).
tsp_order_output_file_path = 'tsp_order.txt'  # Use actual file path

tsp_to_file(dm_file_path, is_classical, tsp_order_output_file_path)
```

From a given distance matrix saved in a file format that supports `numpy.loadtxt`,
you can find its TSP solution either by D-Wave quantum annealer:
```python
is_classical = False
```
or classical algorithm ([python-tsp](https://github.com/fillipe-gsm/python-tsp)):
```python
is_classical = True
```


# Step 3. Construct a MEMC guide tree

`memc_guide_tree.py` takes three arguments: 
(1) a path to a distance matrix, (2) a path to a TSP solution, and (3) a path to save the output.

```python
from memc_guide_tree import MEMCGuideTree

dm_path = 'dm.txt'  # Use actual file path
tsp_path = 'tsp.txt'  # Use actual file path
output_path = 'output.dnd'  # Use actual file path

memc = MEMCGuideTree(dm_path, tsp_path, output_path)
memc.guide_tree_to_file()  # Save the guide tree.
```

memc_guide_tree requires `unionfind.py` implemented by [Debajyoti Nandi](https://github.com/deehzee/unionfind).
`unionfind.py` is under the MIT license and its copy is provided in this repository with the same file name.


---
If you want to build a classical guide tree to compare with our tree, you can use `classic_guide_tree.py`.
```python
from classic_guide_tree import GuideTree

algorithm = "upgma"  # Choose a classical algorithm.
dm_path = "dm.txt"  # Use actual file path
output_path = "output.dnd"  # Use actual file path

gt = GuideTree(dm_path, algorithm, output_path)
gt.guide_tree_to_file()  # Save the guide tree.
```

This code provides the following classical algorithms for guide trees.
All these implementations are based on Biopython with modifications by the author.
- Complete linkage [1]
- Flexible (beta = 0.5) [1]
- Modified UPGMA (in MAFFT) [2]
- Single linkage [1]
- UPGMA [1]
- UPGMC [1]
- Ward [1]
- WPGMA [1]
- WPGMC [1]

References to classical algorithms: 

[1] [Milligan, G.W. Ultrametric hierarchical clustering algorithms. Psychometrika 44, 343-346 (1979)](https://doi.org/10.1007/BF02294699)

[2] [Kazutaka Katoh, John Rozewicki, Kazunori D Yamada, MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization, Briefings in Bioinformatics, Volume 20, Issue 4, July 2019, Pages 1160-1166](https://doi.org/10.1093/bib/bbx108)



# Dependencies
This project requires the following packages:
- [numpy](https://github.com/numpy/numpy)
- [networkx](https://github.com/networkx/networkx)
- [dwave-networkx](https://github.com/dwavesystems/dwave-networkx)
- [dwave-system](https://github.com/dwavesystems/dwave-system)
- [BioPython](https://github.com/biopython/biopython)
- [python-tsp](https://github.com/fillipe-gsm/python-tsp)

We appreciate contributors and maintainers of these valuable packages.



# License
Copyright (c) 2024, Youngjun Park (pyoungj09 at gmail dot com)

Each code file includes a header explaining the applicable license.
Please see the LICENSE file in the project root for more information.



# Citing this work
(Link will be added).