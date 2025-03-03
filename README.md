# compare_trees_RF
compare phylogenetic trees and calculate TMRCAs of recombinant forms 


## get_RF_halflife.py

Usage
```
get_RF_halflife.py [-h] -tree1 TREE_BEAST -tree2 TREE2 [-pthr POSTERIOR_THRESHOLD]

options:
  -h, --help            show this help message and exit
  -tree1 TREE_BEAST, --tree_beast TREE_BEAST
                        Path to time tree in nexus format inferred using BEAST software
  -tree2 TREE2, --tree2 TREE2
                        Path to phylogenetic tree in nwk or nexus format
  -pthr POSTERIOR_THRESHOLD, --posterior_threshold POSTERIOR_THRESHOLD
                        Threshold for posterior values of nodes that to count. Ranges from 0 to 1.
```


## get_subtrees.R

## plot_trees_test.R
