library(ggtree)
library(treeio) #read.beast
library(ape) # for phylo objects
library(tidyr)
library(ggpubr)

# set  working directory
#setwd("")
# Load functions from source code
source("get_subtrees.R")

# Example 1. The nexus trees, leaf sets are different. 
#The set of leaves in one tree is a subset of the leaves in another tree

# converts file with common subtrees in newick format into list with taxa names
# *_commontrees.txt" is a file with common subtrees between two phylogenetic trees
# Each subtree is in newick format
# this file is output of get_RF_halflife.py
# Example: python.exe .\get_RF_halflife.py -tree1 .\norovirus_vp1.tree -tree2 .\norovirus_rdrp_g2.tree

subtrees = get_subtrees("norovirus_vp1_norovirus_rdrp_g2_commontrees.txt")

# Path to file with tree in nexus format produced by BEAST
tree1_name = "norovirus_vp1.tree"
tree1 = read.beast(tree1_name)

tree2_name = "norovirus_rdrp_g2.tree"
tree2 = read.beast(tree2_name)

# Plot beast subtree and color braches of common subtrees
p1 = ggtree(tree1,mrsd="2025-01-01") +  geom_tiplab(size=2) + 
            geom_label2(aes(label=round(as.numeric(posterior),2)),size=2)
p1 = groupOTU(p1, subtrees, 'Clade') + aes(color=Clade) +
  theme(legend.position="right") + scale_color_manual(values=c("black", "firebrick"))

p2 = ggtree(tree2,mrsd="2025-01-01") +  geom_tiplab(size=2) + 
            geom_label2(aes(label=round(as.numeric(posterior),2)),size=2)
p2 = groupOTU(p2, subtrees, 'Clade') + aes(color=Clade) +
  theme(legend.position="right") + scale_color_manual(values=c("black", "firebrick"))

ggarrange(p1,p2)

# Example 2. Trees with the same leaf sets


# Path to file with tree in nexus format produced by BEAST
tree1_name = "Asia1_P1.tree"
tree1 = read.beast(tree1_name)

tree2_name = "Asia1_Coord_4981_6100.treefile"
tree2 = read.tree(tree2_name)

# common subtrees determined from coinciding bipartitions in two trees
subtrees = get_subtrees("Asia1_P1_Asia1_Coord_4981_6100_commontrees_bip.txt")


# Plot beast subtree and color braches of common subtrees
p1 = ggtree(tree1,mrsd="2025-01-01") +  geom_tiplab(size=2) + 
      geom_label2(aes(label=round(as.numeric(posterior),2)),size=2) +
      theme_tree2()
p1 = groupOTU(p1, subtrees, 'Clade') + aes(color=Clade) +
  theme(legend.position="right") + scale_color_manual(values=c("black", "firebrick")) 

p2 = ggtree(tree2) +  geom_tiplab(size=2) + 
      geom_label2(aes(label=round(as.numeric(label),2)),size=2) + 
      theme_tree2()
p2 = groupOTU(p2, subtrees, 'Clade') + aes(color=Clade) +
  theme(legend.position="right") + scale_color_manual(values=c("black", "firebrick"))

ggarrange(p1,p2)
