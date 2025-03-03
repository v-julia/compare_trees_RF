library(ggtree)
library(treeio) #read.beast
library(ape) # for phylo objects
library(tidyr)
source("get_subtrees.R")

subtrees = get_subtrees("trees.txt")

tree1_name = "tree1.tree"
tree1_name = "tree2.tree"
tree1 = read.beast(tree1_name)

p <- ggtree(tree1,mrsd="2025-01-01")
groupOTU(p, subtrees, 'Clade') + aes(color=Clade) +
  theme(legend.position="right") + scale_color_manual(values=c("black", "firebrick"))

