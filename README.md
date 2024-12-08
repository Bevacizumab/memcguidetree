# The Minimum Evolution with the Molecular Clock (MEMC) Guide Tree Algorithm
This repository contains the source code for the distance-based guide tree construction algorithm introduced in [our paper](https://doi.org/10.1101/2024.11.30.626202).
Distance-based guide trees are constructed using [hierarchical clustering algorithms](https://en.wikipedia.org/wiki/Hierarchical_clustering) and are utilized in [multiple sequence alignments (MSAs)](https://en.wikipedia.org/wiki/Multiple_sequence_alignment).
These trees determine the order in which input sequences are aligned.

The minimum evolution (ME) is associated with minimizing the total sum of branch lengths. [The problem of finding the optimal tree corresponding to the ME can be mapped to a traveling salesman problem (TSP).](https://doi.org/10.1093/bioinformatics/16.7.619)
In this project, quantum annealing is utilized to solve these TSPs with scalability.
Our algorithm (`src/memc_guide_tree.py`) constructs a guide tree from a given TSP solution with linear time complexity (*O(N)*).
The molecular clock hypothesis is employed in this construction process.

For simplicity, we have adopted the sequence labeling convention used in [MAFFT](https://mafft.cbrc.jp/alignment/software/); sequences in a guide tree are labeled as 1, 2, ... *n*.
To input this guide tree into MAFFT, it must be converted to MAFFT's native tree format using a [Ruby script](https://mafft.cbrc.jp/alignment/software/newick2mafft.rb).

This code also applies to various MSA tools that utilize distance-based guide trees, with simple modifications to the sequence labeling convention.
We encourage you to make the necessary modifications and use this code for other MSA tools as long as you adhere to the included LICENSE file and provide proper citation of the code and [our paper](https://doi.org/10.1101/2024.11.30.626202).

Tests for this guide tree algorithm can be performed using various MSA benchmarks.  
To use [BAliBASE 3.0](https://www.lbgi.fr/balibase/), you can calculate scores using [this Git repository](https://github.com/robinhundt/bali-score) (a Rust compiler is required here, which can be installed via [rustup](https://rustup.rs/)).




# Example Usage:

This project provides an implementation for constructing guide trees compatible with [MAFFT](https://mafft.cbrc.jp/alignment/software/).
As mentioned in the introduction, modifying the sequence labeling convention allows this code to be used for any MSA tool.  
The process involves three main steps:

1. **Obtain a Distance Matrix from an MSA Program**  
2. **Solve a TSP**  
3. **Construct an MEMC Guide Tree**



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
For more details, refer to [the manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html) or [the original paper](https://doi.org/10.1093/molbev/mst010).
Note that `--maxiterate` must be set to 0, and the command should not include an output file directory.  
This command generates a distance matrix in a `.hat2` file, which is created in the same directory as the input file.

A typical output format for four input sequences:
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
Here, `d(i,j)` represents the distance between the *i*-th and the *j*-th sequences.  
MAFFT outputs distances with three decimal points by default.

For more precise guide trees, the `WriteFloatHat2_pointer_halfmtx` function in `io.c` should be modified to print distances with greater precision (at least five or six decimal points).
For instance, replace the following `//Old` lines in the function with `//New` lines and recompile MAFFT to print ten decimal points:
```c
// Old
fprintf( hat2p, " %#6.3f\n", max * 2.5 );

// Old
fprintf( hat2p, DFORMAT, mtx[i][j-i] );
```
```c
// New
fprintf( hat2p, " %#6.10f\n", max * 2.5 );

// New
fprintf( hat2p, " %6.10f", mtx[i][j-i] );
```

`memc_guide_tree.py` and `classic_guide_tree.py` scripts expect a distance matrix saved in the following format:
```
0.00000 d(1,2)  d(1,3)  d(1,4)
d(1,2)  0.00000 d(2,3)  d(2,4)
d(1,3)  d(2,3)  0.00000 d(3,4)
d(1,4)  d(2,4)  d(3,4)  0.00000
```
We do not provide the code for extracting the distance information from `.hat2` files and converting it to the above format.
This is because the process is straightforward, and distance matrix outputs cannot be generalized across various MSA tools.


## Step 2. Solve a TSP

Use `tsp_from_dm.py` to solve the TSP:
```python
from tsp_from_dm.py import *

dm_file_path = 'dm.txt'  # Specify the actual file path
is_classical = False  # D-Wave Quantum Annealer (LeapHybridSampler).
tsp_order_output_file_path = 'tsp_order.txt'  # Specify the actual file path

tsp_to_file(dm_file_path, is_classical, tsp_order_output_file_path)
```
[D-Wave quantum annealer](https://www.dwavesys.com/), especially `LeapHybridSampler` is used in this project.
Given a distance matrix saved in a file format that supports `numpy.loadtxt`,
you can solve the TSP using either D-Wave quantum annealer:
```python
is_classical = False
```
or a classical algorithm ([python-tsp](https://github.com/fillipe-gsm/python-tsp)):
```python
is_classical = True
```


## Step 3. Construct an MEMC guide tree

`memc_guide_tree.py` requires three arguments: 
(1) a path to the distance matrix, (2) a path to the TSP solution, and (3) a path to save the output.

```python
from memc_guide_tree import MEMCGuideTree

dm_path = 'dm.txt'  # Specify the actual file path
tsp_path = 'tsp.txt'  # Specify the actual file path
output_path = 'output.dnd'  # Specify the actual file path

memc = MEMCGuideTree(dm_path, tsp_path, output_path)
memc.guide_tree_to_file()  # Save the guide tree.
```

`memc_guide_tree.py` requires `unionfind.py`, implemented by [Debajyoti Nandi](https://github.com/deehzee/unionfind).
`unionfind.py` is distributed under the MIT license, and a copy is included in this repository with the same file name.

If you want to use this code to generate guide trees for other MSA tools, modify the `_data_load` method in the `MEMCGuideTree` class to load the appropriate sequence labels.
However, proceed with caution, as modifications may introduce errors.
For example, `UnionFind` may not function properly in this code if sequence labels are not integers, as this scenario has not been tested.


---
If you want to build a classical guide tree to compare with our tree, you can use `classic_guide_tree.py`:
```python
from classic_guide_tree import GuideTree

algorithm = "upgma"  # Choose a classical algorithm.
dm_path = "dm.txt"  # Specify the actual file path
output_path = "output.dnd"  # Specify the actual file path

gt = GuideTree(dm_path, algorithm, output_path)
gt.guide_tree_to_file()  # Save the guide tree.
```

This code provides the following classical algorithms for constructing guide trees.
All these implementations are based on Biopython with modifications by the author.
- Complete linkage [1]
- Flexible (beta = 0.5) [1]
- Modified UPGMA (default in MAFFT) [2]
- Single linkage [1]
- UPGMA [1]
- UPGMC [1]
- Ward [1]
- WPGMA [1]
- WPGMC [1]

To use a specific algorithm, set the `algorithm` variable to one of the strings in the following list:
```
['complete_linkage', 'flexible', 'modified_upgma', 'single_linkage', 'upgma', 'upgmc', 'ward', 'wpgma', 'wpgmc']
```

References to the classical algorithms: 

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

We greatly appreciate the contributors and maintainers of these valuable packages.
Specific versions of these packages are not specified.



# License
Copyright (c) 2024, Youngjun Park (pyoungj09 at gmail dot com)

Each code file contains a header explaining the applicable license.
For more details, please refer to the LICENSE file in the project root.


# Citing this work
If you find this project useful, please cite our paper:

```bibtex
@article {Park2024.11.30.626202,
	author = {Park, Youngjun and Kim, Juhyeon and Huh, Joonsuk},
	title = {Scalable Guide Tree Construction Using Quantum Annealing for Multiple Sequence Alignment},
	year = {2024},
	doi = {10.1101/2024.11.30.626202},
	publisher = {Cold Spring Harbor Laboratory},
	journal = {bioRxiv}
}
```