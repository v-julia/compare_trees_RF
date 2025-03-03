library(ggtree)
library(treeio) #read.beast
library(ape) # for phylo objects
library(tidyr)
source("get_subtrees.R")

setwd("D:/Virology/RNA_viruses/Compare_trees")

subtrees = get_subtrees("trees.txt")

tree1_name = "D:/Virology/RNA_viruses/Picornaviruses/virus-recombination/Trees/MCC/Asia1_P1.tree"
tree1_name = "norovirus_rdrp_g2_yearonly_colnodes.tree"
tree1 = read.beast(tree1_name)

p <- ggtree(tree1,mrsd="2025-01-01")
groupOTU(p, subtrees, 'Clade') + aes(color=Clade) +
  theme(legend.position="right") + scale_color_manual(values=c("black", "firebrick"))

