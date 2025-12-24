# RF-HL (calculate Recombinant Forms' Half-life)

RF-HL is a tool for estimating the half-life of virus recombinant forms. 

This repository contains functions for identifying and visualizing common subtrees in two phylogenetic trees. One of the trees should be a time-scaled phylogeny produced by BEAST.


## get_RF_halflife.py
This script provides two methods to compare two trees and calculate median ages of coinciding nodes ("method" parameter).

1) subtrees
Finds common subtrees between two trees and calculates median of their MRCA ages. The nodes' tMRCA are extracted from the first tree. This tree should be a time-scaled phylogeny produced by BEAST software. The second tree's format could be nexus or newick. Output:
- *_commontrees.txt - text file with common subtrees in newick format
- *_heights.txt - ages of common subtrees
Importantly, this method allows the trees' taxa not to coincide. 


<img src="https://github.com/v-julia/compare_trees_RF/blob/main/images/Common_subtrees_scheme.png" align="center" width=1200/>

Figure 1. Scheme of "subtrees" method.

2) bipartitions
Finds common bipartitions, or non-trivial branches between two trees, calculates median ages of nodes that correspond to bipartitions. The nodes' tMRCA are extracted from the first tree. This tree should be a time-scaled phylogeny produced by BEAST software. The second tree's format could be nexus or newick. Output:
- *_commontrees_bip.txt - text file with common subtrees in newick format. Nested nodes are removed.
- *_heights_bip.txt - ages of common subtrees. Nested nodes are removed.
Taxa in two trees must coincide. 

<img src="https://github.com/v-julia/compare_trees_RF/blob/main/images/Common_bipartitions_scheme.png" align="center" width=1200/>

Figure 2. Scheme of "bipartitions" method.

3) all
Applies both methods. Output a table with the following columns:
```
- tree1 - Name of file with the first tree without file extension
- tree2 - Name of file with the second tree without file extension
- bipartitions - number of non-trivial branches in the first tree
- bipartitions1 posterior>threshold - number of non-trivial branches in the first tree with posterior higher than a threshold
- bipartitions2 bootstrap>threshold2 - number of non-trivial branches in the second tree with posterior/bootstrap support higher than a threshold
- coinciding bipartitions no threshold - number of coinciding non-trivial branches in two trees
- coinciding bipartitions - number of coinciding non-trivial branches with posterior/bootstrap support higher than threshold in two trees
- RF times/bipartitions no threshold - median ages of recombinant forms (common subtrees) calculated from common bipartitions 
- RF times/bipartitions - median ages of recombinant forms (common subtrees) calculated from common bipartitions with posterior/bootstrap support higher than a threshold
- coinciding bipartitions no threshold (no nested clades) - number of coinciding non-trivial branches in two trees, nested clades are deleted
- coinciding bipartitions (no nested) -  number of coinciding non-trivial branches with posterior/bootstrap support higher than threshold in two trees, nested clades are deleted
- RF times/bipartitions no threshold (no nested) - median ages of recombinant forms (common subtrees) calculated from common bipartitions, nested clades are deleted
- RF times/bipartitions (no nested) - median ages of recombinant forms (common subtrees) calculated from common bipartitions with posterior/bootstrap support higher than threshold, nested clades are deleted
- coinciding subtrees no threshold - number of coinciding subtrees
- RF times/subtrees no threshold - median ages of recombinant forms (common subtrees) 
- coinciding subtrees - number of coinciding subtrees with posterior supports in the first tree higher than a threshold
- RF times/subtrees - median ages of recombinant forms (common subtrees) with posterior supports in the first tree higher than a threshold
- coinciding subtrees (thr while encoding) - number of coinciding subtrees with posterior supports in the first tree higher than a threshold. The nodes with supports lower than `thr` and `thr2` where not considered while encoding the tree.
- RF times/subtrees(thr while encoding) - median ages of recombinant forms (common subtrees) with posterior supports in the first tree higher than a threshold. The nodes with supports lower than `thr` and `thr2` where not considered while encoding the tree.


```


Usage
```
get_RF_halflife.py [-h] -tree1 TREE_BEAST -tree2 TREE2 -method METHOD [-thr POSTERIOR_THRESHOLD]
                          [-thr2 THRESHOLD2] [-out OUTPUT_DIR]

options:
  -h, --help            show this help message and exit
  -tree1 TREE_BEAST, --tree_beast TREE_BEAST
                        Path to time tree in nexus format inferred using BEAST software
  -tree2 TREE2, --tree2 TREE2
                        Path to phylogenetic tree in nwk or nexus format
  -method METHOD, --method METHOD
                        Method for determining recombinant forms. 'subtrees', 'bipartitions' or 'all'. If method==all,
                        outputs a table with numbers of coinciding partitions and subtrees, and ages of RFs half-lives
  -thr POSTERIOR_THRESHOLD, --posterior_threshold POSTERIOR_THRESHOLD
                        Threshold for posterior values of nodes to count. Ranges from 0 to 1.
  -thr2 THRESHOLD2, --threshold2 THRESHOLD2
                        Threshold for bootstrap values or posterior values of branches to count (for 'bipartitions'
                        method). Bootstrap values range from 0 to 100. Posterior probabilities range from 0 to 1
  -out OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory
```


## get_subtrees.R

Helper functions for parsing *_commontrees.txt


## plot_trees_test.R

Example of visualization.
