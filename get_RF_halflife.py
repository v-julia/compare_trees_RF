import argparse
import hashlib
import re
import sys
import os

import dendropy as dpy
import numpy as np
import pandas as pd

# Regular expression for decimals
locval_re = re.compile(r'[0-9]+\.[0-9E\-]+')

'''
This function traverses beast tree and adds the attributes 'posterior' and 'height' for all tree nodes according to the annotations.
Parameters
----------
tree: class dendropy.datamodel.treemodel.Tree object
    A tree inferred by BEAST

Returns
-------
tree: class dendropy.datamodel.treemodel.Tree object
    Tree with updated nodes with 'posterior' and 'height' attributes.
'''

def parse_beast_tree_node_info(tree):

    for nd in tree.postorder_node_iter():
        node_comment = nd._get_annotations()
        for annot in node_comment:
            if str(annot).startswith('posterior='):
                posterior =  float(locval_re.search(str(annot)).group())
                nd.posterior = posterior
            if str(annot).startswith('height='):
                height =  float(locval_re.search(str(annot)).group())
                nd.height = height

    return tree


'''
This function sorts alphabetically taxa names in bipartition produced by split_as_newick_string() (dendropy)

Parameters
----------
bipart_str: string
    Bipartition split represented as a newick string

Returns
-------
bipart_str: string
    bipart_str with sorted taxa

Example
-------
bipart_str = "((B,A),(C));"

sort_bipart(bipart_str)
"((A,B),(C));"
'''
def sort_bipart(bipart_str):
    if '((' in bipart_str:
        bipart_l = bipart_str.split('), (')
        #print(bipart_l)
        for i in range(len(bipart_l)):
            bipart_l[i] = bipart_l[i].strip('((').strip('));')
            if len(bipart_l[i]) == 1:
                continue
            else:
                part_l = bipart_l[i].split(', ')
                part_l.sort()
                bipart_l[i] = ','.join(part_l)
        
        #bipart_l[0] = bipart_l[0].strip('((')
        #bipart_l[1] = bipart_l[1].strip('))')
    
        #print(bipart_l)
        bipart_l.sort()
        #print(bipart_l)
    
    
        bipart_str = '((' + bipart_l[0] + ')' + ',' + '(' + bipart_l[1] + '))'
        return bipart_str
    else:
        bipart_sorted = bipart_str.strip("(").strip(");").split(",")
        bipart_sorted.sort()
        bipart_str = "(" + ",".join(bipart_sorted) +  ")"
        #return ""
        return bipart_str


'''
Calculates hash digest (sha256) for a string
Parameters
----------
stri: string

Returns
-------
stri: string
    digest (sha256) of stri as a string object of double length, containing only hexadecimal digits. 
''' 
def stri2hash(stri):
    h = hashlib.new('sha256',usedforsecurity=False)
    #h.update(json.dumps(nd._as_newick_string(edge_lengths=None)).encode('utf-8'))
    h.update((stri).encode('utf-8'))
    return h.hexdigest()

'''
A simple function for calculating a hash for a node. 
Represents a tree node as a string with enumeration of all possible bipartitions of node's subtree. Returns a hash of such a string and the string itself.
This function does not consider the support of subtree nodes.

Parameters
---------------
node - class dendropy.datamodel.treemodel.Node
    Node of a tree

Output
---------------
stri_bipart: string 
    Node representation as a string with enumeration of all possible bipartitions of node's subtree

stri2hash(stri_bipart): string
    Hash digest for stri_bipart as a string object of double length, containing only hexadecimal digits.
'''
def node2hash(node):
    subtree = dpy.Tree(seed_node=node.extract_subtree())
    biparts = subtree.encode_bipartitions(collapse_unrooted_basal_bifurcation=False)

    bipart_list = [sort_bipart(bip.split_as_newick_string(subtree.taxon_namespace)) for bip in biparts]
    bipart_list.sort()
    stri_bipart = ';'.join(bipart_list[1:])

    return stri2hash(stri_bipart),stri_bipart

'''
Represents a tree as a dictionary where 
keys are nodes encoded as hash digest for a string with enumeration of all possible bipartitions of node's subtree,
values are supports (bootstrap, posterior probability etc) of the nodes

Parameters
---------------
tree: class dendropy.datamodel.treemodel.Tree object

treetype: string(default='beast')
    Type of a tree. If 'beast', the tree is supposed to be inferred by beast, nodes of such tree are supposed to have 'posterior' attribute.


Returns
-------
node2support: dict
    A dictionary with support values for inner nodes of the tree.
    keys - nodes encoded as hash digest for a string with enumeration of all possible bipartitions of node's subtree,
    values - supports (bootstrap, posterior probability etc) of the nodes
    node2support[node hash] = support

'''
def create_node2support_dict(tree, treetype="beast"):
    node2support = {}
    
    for nd in tree.postorder_node_iter():
        #if nd.is_leaf():
        #    continue
        if treetype != "beast":
            if nd.label!= None:
                support = float(nd.label)
                
            else:
                # Case for non-binary tree, when several leafs cpme from one node
                support = 0
                #print(nd)
        else:
            if hasattr(nd, 'posterior'):
                support = nd.posterior
            else:
                support = 0
                #print(nd)
        #print(node2hash(nd))
        nd_hash, stri = node2hash(nd)
        node2support[node2hash(nd)] = support

    return node2support

'''
A complicated function for calculating a hash for a node. 
Represents a tree node as a string with enumeration of bipartitions of node's subtree with support higher than user-defined threshold. 
Returns a hash of such a string and the string itself.

Parameters
---------------
node: class dendropy.datamodel.treemodel.Node
    Node of a tree
node2support: dict 
    A dictionary with support values for inner nodes of the tree.
    keys - nodes encoded as hash digest for a string with enumeration of all possible bipartitions of node's subtree,
    values - supports (bootstrap, posterior probability etc) of the nodes
thr: float
    Threshold for node support
Output
---------------
stri_bipart: string 
    Node representation as a string with enumeration of all possible bipartitions of node's subtree

stri2hash(stri_bipart): string
    Hash digest for stri_bipart as a string object of double length, containing only hexadecimal digits.
'''

def node2hash_thr(node, node2support, thr):
    subtree = dpy.Tree(seed_node=node.extract_subtree())
    subtree.encode_bipartitions(collapse_unrooted_basal_bifurcation=False)
    bipart_list_thr = []
    #print("New node")
    #print(node)
    #print("Node nwk")
    #print(subtree._as_newick_string(edge_lengths=None))
    #print("Namespace")
    node_labels = subtree.taxon_namespace.labels()
    node_labels.sort()
    #print(subtree.taxon_namespace.labels())
    for edge in subtree.levelorder_edge_iter():
        if edge.is_leaf():
            continue
        #print(edge.bipartition.split_as_newick_string(subtree.taxon_namespace))
        bipart_sorted = sort_bipart(edge.bipartition.split_as_newick_string(subtree.taxon_namespace)).strip(";")
        #print("Bipartition")
        #print(bipart_sorted)

        #print("Node2hash")
        #print(node2hash(edge.head_node))
        support = node2support[node2hash(edge.head_node)]
        #print(support)
        if support > thr:
            #print("Added")
            bipart_list_thr.append(bipart_sorted)
            #print("Edge nwk")
            #print(bipart.split_as_newick_string(subtree.taxon_namespace))
    #print(bipart_list_thr)
    if len(bipart_list_thr) == 0:
        stri_bipart = "(" + ",".join(node_labels) + ")*"
    else:
        stri_bipart = ';'.join(bipart_list_thr)
    #print("Final str bipart")
    #print(stri_bipart)
    return stri2hash(stri_bipart), stri_bipart

'''
Postorder traverses the tree and adds the attribute 'hash' for each node. 
Hash (sha256) is calculated for string that with enumeration of bipartitions with support higher than user-defined threshold in subtree with the descendants of a particular node.
Returns dictionary with hashes and corresponding subtrees in newick format.


Parameters
----------
tree: class dendropy.datamodel.treemodel.Tree
Phylogenetic tree

thr: float
    Threshold support values in inner nodes of subtree with the descendants of a particular node. 
    If the support of inner node is lower than threshold, the corresponding bipartition will not be added to a the list of bipartitions that represent the node.

timescaled: bool, default=True
    If True, the tree is treated as a time-tree inferred by BEAST. For each node, the height of the youngest leaf in the corresponding subtree is stored.

treetype: string (default="beast")
    Type of input tree. If "beast", the time tree inferred by BEAST is supposed.

Returns
-------
tree_hashes: dict
Dictionary:
key - node encoded as hash digest for a string with enumeration of all possible bipartitions of node's subtree (inner node supports are not considered)
value - newick substring for subtree that corresponds to the node

subtree_heights: dict
key - node encoded as hash digest for a string with enumeration of all possible bipartitions of node's subtree (inner node supports are not considered)
value - height of the youngest leaf of subtree that corresponds to the node

#biparts: dict
#Dictionary with corresponding of node hash and string with enumeration of all possible bipartitions of node's subtree
'''

def add_hashes(tree, thr=None, timescaled=False, treetype="beast"):
    tree_hashes = {}
    #biparts = {}
    nd2support = create_node2support_dict(tree, treetype)
    
    if timescaled:
        subtree_heights = {}

    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            continue
        else:
            if thr != None:
                nd.hash, stri_bipart = node2hash_thr(nd, nd2support, thr)
                #print(stri_bipart)
            else:
                nd.hash, stri_bipart = node2hash(nd)
            
            tree_hashes[nd.hash] = nd._as_newick_string(edge_lengths=None)
            tree_hashes[nd.hash] = re.sub("\)[0-9]+", ")", tree_hashes[nd.hash])
            #biparts[nd.hash] = stri_bipart
            
            if timescaled:
                min_leaf_height = 100
                for leaf in nd.leaf_iter():
                    if leaf.height < min_leaf_height:
                        min_leaf_height = leaf.height
                subtree_heights[nd.hash] = min_leaf_height

    if timescaled:
        return tree_hashes, subtree_heights#, biparts
    return tree_hashes#, biparts

'''
Extracts heights of nodes which subtrees correspond in two trees

Parameters
----------
tree1: class dendropy.datamodel.treemodel.Tree object
    Time inferred using BEAST. The nodes of this tree are supposed to have attribute height and hash. 
    The attribute posterior is also needed if the threshold for posterior value is used.

subtree_times1: dict
    Dictionary with hashes and heights of the youngest tip in tree1

hashes_tree2: dict
    Dictionary with hashes for the second tree you would compare with

posterior_thr: float
    Threshold posterior value. If posterior of the node is lower that the threshold, the height of subtree will not be counted.

Returns
-------
heights: list
    List with heights of subtrees common between tree1 and tree2. Height is calculated as difference between the age of the youngest leaf and node height.

'''
def get_common_subtrees(tree1, subtree_times1, hashes_tree2, posterior_thr = None):
    heights = []
    common_subtrees = []
    stack = [tree1.seed_node]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            continue
        if node.hash in hashes_tree2.keys():
            # check node posterior
            if posterior_thr:
                if node.posterior > posterior_thr:
                    #yield node.height
                    #print(node.height)
                    heights.append(node.height - subtree_times1[node.hash])
                    common_subtrees.append(hashes_tree2[node.hash])
                else:
                    stack.extend(n for n in reversed(node._child_nodes))
            else:
                heights.append(node.height - subtree_times1[node.hash])
                common_subtrees.append(hashes_tree2[node.hash])
        else:
            stack.extend(n for n in reversed(node._child_nodes))
    return heights, common_subtrees




'''
Splits tree to bipartitions. For newick files, saves bootstrap support of a corresponding bipartitions. For nexus files obtained from BEAST saves posterior support and node height.

Parameters
----------
tree: class dendropy.datamodel.treemodel.Tree object
    Tree in newick format or in nexus format (inferred using BEAST).

treetype: string (default='newick')
    Type of tree ('newick' or 'nexus')

Returns
-------
bipartitions: dict
    A dictionary of tree bipartitions where:
        keys - hash digest (as string) for a bipartition encoded as newick string
        values - list.
    For newick tree, bipartitions[key]=[bootstrap, sorted bipartition as newick string]
    For nexus tree bipartitions[key]=[posterior, node height, node height - min_leaf_height, sorted bipartition, newick string for a corresponding node]

'''

def encode_bipartitions(tree, treetype="newick"):
    #counter for non binary nodes
    k=0
    bipartitions = {}
    tree.encode_bipartitions(suppress_storage=True)
    for bipart in tree.bipartition_edge_map:
        edge = tree.bipartition_edge_map[bipart]
        #print(edge.head_node)
        #print(bipart.split_as_newick_string(tree.taxon_namespace))
        if edge.is_leaf():
            continue
        else:
            bipart_sorted = sort_bipart(bipart.split_as_newick_string(tree.taxon_namespace))
            h = hashlib.new('sha256',usedforsecurity=False)
            h.update((bipart_sorted).encode('utf-8'))
            if treetype == "newick":
                if edge.head_node.label!= None:
                    support = float(edge.head_node.label)
                else:
                    # Case for non-binary tree, when several leafs cpme from one node
                    support = 0
                    k+=1
    
                #print(bipart.split_as_newick_string(tree.taxon_namespace))
                bipartitions[h.hexdigest()] = [support,bipart_sorted]
            else:
                support = edge.head_node.posterior
                height = edge.head_node.height

                min_leaf_height = 100
                for leaf in edge.head_node.leaf_iter():
                    if leaf.height < min_leaf_height:
                        min_leaf_height = leaf.height
                newick_substr = re.sub("\)[0-9]+", ")", edge.head_node._as_newick_string(edge_lengths=None))  
                bipartitions[h.hexdigest()] = [support, height, height - min_leaf_height,bipart_sorted,newick_substr]
    print("Bipartitions from non-binary nodes = {}".format(k))
    return(bipartitions)

'''
Extracts heights of nodes which bipartitions correspond in two trees

Parameters
----------
biparts_tree1: dict
Dictionary for bipartitions of tree1 which is time tree in nexus format inferred using BEAST. 
    biparts_tree1[hash] = [support, height, height - min_leaf_height, bipart_sorted, newick_substr], where 
            support - is a support value of node corresponding to a bipartition,
            height is a node heigth, bipart_sorted is a string with sorted bpartitions,
            newick_substr is a newick substring for a corresponding node

biparts_tree2: dict
    Dictionary for bipartitions of tree2 which can be both in newick and time tree in nexus format

Taxa labels must coincide in two trees!
The times and common subtrees are calculated using tree1.

Returns
-------
heights: list
    List with heights of nodes that correspond to bipartitions common in tree1 and tree2

subtrees: list 
    List with subtrees of common bipartitions in tree1 and tree2
'''

def get_common_biparts_old(biparts_tree1, biparts_tree2, posterior_thr=None, bootstrap_sup=None):
    heights = []
    subtrees = []
    for h in biparts_tree1:
        if posterior_thr!=None and biparts_tree1[h][0] < posterior_thr:
            continue
        if h in biparts_tree2.keys():
            #print(biparts_tree1[h][0])
            if bootstrap_sup!=None and biparts_tree2[h][0] < bootstrap_sup:
                continue
            #print(biparts_tree1[h][1])
            #print(biparts_tree1[h][2])
            #print(biparts_tree2[h][0])
            heights.append(biparts_tree1[h][2])
            subtrees.append(biparts_tree1[h][4])

    return heights, subtrees

'''
Extracts heights of nodes which bipartitions correspond in two trees. Excludes nested bipartitions.

Parameters
----------
biparts_tree1: dict
    Dictionary for bipartitions of tree1 which is time tree in nexus format inferred using BEAST. 
    biparts_tree1[hash] = [support, height, height - min_leaf_height, bipart_sorted, newick_substr], where 
            support - is a support value of node corresponding to a bipartition,
            height is a node heigth, bipart_sorted is a string with sorted bpartitions,
            newick_substr is a newick substring for a corresponding node

biparts_tree2: dict
    Dictionary for bipartitions of tree2 which can be both in newick and time tree in nexus format

Taxa labels must coincide in two trees!
The times and common subtrees are calculated using tree1.

Returns
-------
heights: list
    List with heights of nodes that correspond to bipartitions common in tree1 and tree2

subtrees: list 
    List with subtrees of common bipartitions in tree1 and tree2

'''

def get_common_biparts(biparts_tree1, biparts_tree2, posterior_thr=None, bootstrap_sup=None):
    # bipartitions that coincide in two trees
    biparts_coinc = {}
    # subtrees for bipartitions that coincide in two trees
    subtree2hash = {}

    # compare bipartitions
    for h in biparts_tree1:
        if posterior_thr!=None and biparts_tree1[h][0] < posterior_thr:
            continue

        if h in biparts_tree2.keys():
            if bootstrap_sup!=None and biparts_tree2[h][0] < bootstrap_sup:
                continue
            biparts_coinc[h] = biparts_tree1[h]
            subtree2hash[biparts_tree1[h][4]] = h

    # some bipartitions correspond to nested nodes
    # does not consider ages and subtrees for splits of nodes that are descendants of other nodes matching in two trees
    heights = []
    subtrees = []
    for h in biparts_coinc:
        if len(list(filter(lambda x: biparts_coinc[h][4] in x, subtree2hash.keys()))) < 2:
            heights.append(biparts_coinc[h][2])
            subtrees.append(biparts_coinc[h][4])

    return heights, subtrees


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-tree1", "--tree_beast", type=str,
                        help="Path to time tree in nexus format inferred using BEAST software", required=True)
    parser.add_argument("-tree2", "--tree2", type=str,
                        help="Path to phylogenetic tree in nwk or nexus format", required=True)
    parser.add_argument("-method", "--method", type=str,
                        help="Method for determining recombinant forms. 'subtrees', 'bipartitions' or 'all'.\
                        If method==all, outputs a table with numbers of coinciding partitions and subtrees, and ages of RFs half-lives", required=True)
    parser.add_argument("-thr", "--posterior_threshold", type=float,
                        help="Threshold for posterior values of nodes to count. Ranges from 0 to 1.")
    parser.add_argument("-thr2", "--threshold2", type=float,
                        help="Threshold for bootstrap values or posterior values of branches to count (for 'bipartitions' method). \
                        Bootstrap values range from 0 to 100. Posterior probabilities range from 0 to 1")
    parser.add_argument("-out", "--output_dir", type=str,
                        help="Path to output directory")

    args = parser.parse_args()

    if args.output_dir == None:
        args.output_dir = os.getcwd()
    else:
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

    tree1 = dpy.Tree.get_from_path(args.tree_beast, 'nexus')
    print("Getting annotation of tree nodes...")
    tree1 = parse_beast_tree_node_info(tree1)
    print("Done")
    
    
    try:
        tree2 = dpy.Tree.get_from_path(args.tree2, 'nexus')
        tree2_format = 'nexus'
        tree2 = parse_beast_tree_node_info(tree2)
    except:
        try:
            tree2 = dpy.Tree.get_from_path(args.tree2, 'newick')
            tree2_format = 'newick'
            tree2.reroot_at_midpoint()
        except:
            print("Couldn't read tree2!")
            sys.exit(1)
            
    tree1_name = os.path.splitext(os.path.split(args.tree_beast)[-1])[0]
    tree2_name = os.path.splitext(os.path.split(args.tree2)[-1])[0]
    #print(tree1_name)
    
    if args.method == 'bipartitions' or args.method == 'all':
        print("Encoding bipartitions for tree 1...")
        biparts_tree1 = encode_bipartitions(tree1,"nexus")
        print("Done")
        
        print("Encoding bipartitions for tree 2...")
        if tree2_format == 'newick':
            biparts_tree2 = encode_bipartitions(tree2,"newick")
        else:
            biparts_tree2 = encode_bipartitions(tree2,"nexus")
        print("Done\n")
        
        
        # COMMON BIPARTITIONS (NO THRESHOLD)
        
        biparts_values = []
        # Calculate all common bipartitions and common bipartitions without nested nodes with no threshold for  method check
        for bip_function, text in zip([get_common_biparts_old, get_common_biparts], [""," (no nested nodes)"]):
        
            # get common bipartitions (no threshold)
            heights_bip_nothr, subtrees_bip_nothr = bip_function(biparts_tree1, biparts_tree2)
            
            df = pd.DataFrame(biparts_tree1).T
            df.columns = ['posterior', 'height_raw', 'height_corrected', 'bip', 'subtree']
            if text == "":
                print("Number of bipartitions in tree1 = {}".format(str(df.shape[0])))
                print("Number of bipartitions with posterior > {} in tree1 = {}".format(str(args.posterior_threshold),df[df['posterior']>args.posterior_threshold].shape[0]))

            if tree2_format == 'newick':
                df2 = pd.DataFrame(biparts_tree2).T
                df2.columns = ['support', 'bip']
                bip_thr_tree2 = df2[df2['support']>args.threshold2].shape[0]
                if text == "":
                    print("Number of bipartitions with bootstrap > {} in tree2 = {}".format(str(args.threshold2), bip_thr_tree2))
            else:
                df2 = pd.DataFrame(biparts_tree2).T
                df2.columns = ['posterior', 'height_raw', 'height_corrected', 'bip', 'subtree']
                bip_thr_tree2 = df2[df2['posterior']>args.threshold2].shape[0]
                if text == "":
                    print("Number of bipartitions with posterior > {} in tree2 = {}".format(str(args.threshold2),bip_thr_tree2))
            
            print("\n")
            print("Number of coinciding bipartitions"+ text +" in tree1 and tree2 with no thresholds {}".format(len(subtrees_bip_nothr)))
            print("The median height of common bipartitions"+ text +" with no thresholds is {}".format(round(np.median(heights_bip_nothr),4)))

            print("\n")
            
            # get common bipartitions (with threshold)
            heights_bip, subtrees_bip = bip_function(biparts_tree1, biparts_tree2, args.posterior_threshold, args.threshold2)        
            print("Number of coinciding bipartitions"+ text +" in tree1 and tree2 {}".format(len(subtrees_bip)))
            print("The median height of common bipartitions"+ text +" is {}".format(round(np.median(heights_bip),4)))
            print("\n")
            
            if text == "":
                biparts_values = biparts_values +  [df.shape[0],
                                                    df[df['posterior']>args.posterior_threshold].shape[0],
                                                    bip_thr_tree2,
                                                    len(subtrees_bip_nothr),
                                                    len(subtrees_bip),
                                                    float(round(np.median(heights_bip_nothr),4)),
                                                    float(round(np.median(heights_bip),4))]
            else:
                biparts_values = biparts_values +  [len(subtrees_bip_nothr),
                                                    len(subtrees_bip),
                                                    float(round(np.median(heights_bip_nothr),4)),
                                                    float(round(np.median(heights_bip),4))]
    
            with open(os.path.join(args.output_dir,tree1_name + '_' + tree2_name + "_commontrees_bip"+ text.strip().replace(" ", "-")+".txt"), 'w') as file:
                file.write("\n".join(subtrees_bip) + "\n")
            file.close()

            with open(os.path.join(args.output_dir, tree1_name + '_' + tree2_name + "_heights_bip"+ text.strip().replace(" ", "-")+".txt"), 'w') as file:
                file.write("\n".join([str(h) for h in heights_bip]) + "\n")
            file.close()

    if args.method == 'subtrees' or args.method == 'all':
        # Calculating hashes without considering nodes' supports
        print("Calculating hashes for tree 1...")
        hashes_tree1, subtree_times1 = add_hashes(tree1, timescaled=True)
        print("Done")

        print("Calculating hashes for tree 2...")
        hashes_tree2 = add_hashes(tree2)
        print("Done")

        print("Comparing trees...")
        
        subtrees_values = []
        heights, subtrees = get_common_subtrees(tree1, subtree_times1, hashes_tree2)
        subtrees_values = subtrees_values + [len(subtrees), float(round(np.median(heights),4))]
        
        print("Number of coinciding subtrees (no thresholds) = {}".format(len(subtrees)))
        print("The median height of common subtrees(no thresholds) is {}".format(round(np.median(heights),4)))
        
        
        if args.posterior_threshold:
            # find common subtrees
            heights, subtrees = get_common_subtrees(tree1, subtree_times1, hashes_tree2,float(args.posterior_threshold))
            print("Number of coinciding subtrees = {}".format(len(subtrees)))
            print("The median height of common subtrees is {}".format(round(np.median(heights),4)))
            subtrees_values = subtrees_values + [len(subtrees), float(round(np.median(heights),4))]

        
        with open(os.path.join(args.output_dir,tree1_name + '_' + tree2_name + "_commontrees.txt"), 'w') as file:
            file.write("\n".join(subtrees) + "\n")
        file.close()

        with open(os.path.join(args.output_dir, tree1_name + '_' + tree2_name + "_heights.txt"), 'w') as file:
            file.write("\n".join([str(h) for h in heights]) + "\n")
        file.close()
        
        
        if args.posterior_threshold and args.threshold2:
            print("Calculating hashes for tree 1... (nodes with support < {} are not considered)".format(args.posterior_threshold))
            hashes_tree1, times1 = add_hashes(tree1, thr=args.posterior_threshold, timescaled=True)

            print("Calculating hashes for tree 2... (nodes with support < {} are not considered)".format(args.threshold2))
            hashes_tree2 =  add_hashes(tree2, timescaled=False, thr=args.threshold2, treetype="newick")
            heights, subtrees = get_common_subtrees(tree1, times1, hashes_tree2, args.posterior_threshold)
            
            
            print("Number of coinciding subtrees (node supports were considered when encoding trees) = {}".format(len(subtrees)))
            print("The median height of common subtrees (node supports were considered when encoding trees) is {}".format(round(np.median(heights),4)))
            subtrees_values = subtrees_values + [len(subtrees), float(round(np.median(heights),4))]
            
            with open(os.path.join(args.output_dir,tree1_name + '_' + tree2_name + "_commontrees_thr"+str(args.posterior_threshold)+".txt"), 'w') as file:
                file.write("\n".join(subtrees) + "\n")
            file.close()

            with open(os.path.join(args.output_dir, tree1_name + '_' + tree2_name + "_heights_thr"+str(args.posterior_threshold)+".txt"), 'w') as file:
                file.write("\n".join([str(h) for h in heights]) + "\n")
            file.close()

    if args.method == 'all':
        if tree2_format == 'nexus':
            header = ",".join(["tree1",
                                "tree2",
                                "bipartitions",
                                "bipartitions1 posterior>{}".format(args.posterior_threshold),
                                "bipartitions2 posterior>{}".format(args.posterior_threshold),
                                "coinciding bipartitions no threshold",
                                "coinciding bipartitions",
                                "RF times/bipartitions no threshold",
                                "RF times/bipartitions",
                                "coinciding bipartitions no threshold (no nested clades)",
                                "coinciding bipartitions (no nested)",
                                "RF times/bipartitions no threshold (no nested)",
                                "RF times/bipartitions (no nested)",
                                "coinciding subtrees no threshold",
                                "RF times/subtrees no threshold",
                                "coinciding subtrees",
                                "RF times/subtrees",
                                "coinciding subtrees (thr while encoding)",
                                "RF times/subtrees(thr while encoding)\n"])
        else:
            header = ",".join(["tree1",
                                "tree2",
                                "bipartitions",
                                "bipartitions1 posterior>{}".format(args.posterior_threshold),
                                "bipartitions2 bootstrap>{}".format(args.threshold2),
                                "coinciding bipartitions no threshold",
                                "coinciding bipartitions",
                                "RF bipartitions no threshold",
                                "RF bipartitions",
                                "coinciding bipartitions no threshold (no nested clades)",
                                "coinciding bipartitions (no nested)",
                                "RF bipartitions no threshold (no nested)",
                                "RF bipartitions (no nested)",
                                "coinciding subtrees no threshold",
                                "RF subtrees no threshold",
                                "coinciding subtrees",
                                "RF subtrees",
                                "coinciding subtrees (thr while encoding)",
                                "RF times/subtrees(thr while encoding)\n"])

        list_values = [tree1_name, tree2_name] + biparts_values + subtrees_values

        values = ",".join([str(x) for x in list_values]) + '\n'

        with open(os.path.join(args.output_dir, tree1_name + '_' + tree2_name + "_table.csv"), 'w') as file:
            file.writelines([header, values])