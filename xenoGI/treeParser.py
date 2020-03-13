from copy import deepcopy
#file for parsing the tree into a format tha works for the DTLOR algorithm
# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name (will use end vertex to represent) and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop". 
# Similarly, the dummy edge leading to the root of the host tree must be named
#"hTop"


tree=("root",("left",("HI",(),(),2),("HA",(),(),1),1),('right',(),(),2),3)

def tabulate_names(tree):
    """
    Rename internal nodes from the newick tree 
    """
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name!=None:
            pass
        else:
            clade.name = "Node_"+str(idx)
 
def rerootingPruning(bpTree, locusMap):
    """
    Try assigning syntenic locations to internal nodes 
    in the input (arbitrarily rooted) tree. If a loc can be assigned, we are not
    going to root on branches where both terminals are assigned loc
    """
    locus_map=deepcopy(locusMap)
    postorder=list(bpTree.find_clades(order='postorder'))

    count=0
    for node in postorder:

        if node.name in locus_map: pass

        else:
            child1, child2=[child.name for child in node]
            if locus_map[child1]==locus_map[child2]:
                locus_map[node.name]=locus_map[child1]
               
            else:
                locus_map[node.name]="X"  #undetermined
                count+=1

    print("Percentage of nodes not determined: %.2f"%(count/float(len(postorder))))
    return locus_map


def getTreeDictionary(geneTree,tree_dict):
    """
    get the 4-tup tree format from the dictionary format 
    and a dictionary with key as parent node and value as the children tuple
    """
    root=geneTree[0]
    left=geneTree[1]
    right=geneTree[2]
    if left==():
        return tree_dict  #meaningless length
    else:
        tree_dict[root]=(left[0], right[0])
        tree_dict=getTreeDictionary(left, tree_dict)
        tree_dict=getTreeDictionary(right, tree_dict)
        return tree_dict

def left_subtree(tree):
    return tree[1]
def right_subtree(tree):
    return tree[2]
def is_leaf(tree):
    return tree[1]==()
def get_name(tree):
    return tree[0]

def validRooting(tree, locus_map):
    if locus_map[get_name(tree[1])]==locus_map[get_name(tree[2])] and locus_map[get_name(tree[1])]!="X":  
        #don't root at a branch when both terminals are imputed the same loc
       return False
    else: return True

def get_all_rerootings(tree, locus_map):
    """
    Compute all rerootings of a tree

    Parameters
    --------------------
        tree        -- a four tuple tree
        locus_map   -- a dictionary mapping the gene node name
                        to a syntenic location or "X" if cannot be imputed by parsimony

    Return
    --------------------
        rerootings  -- a list of RLR trees 
    """
    #global memorization variables
    #records all distinctive trees by their left and right subtree roots (branch where rooting happens)
    alltrees=set()
    #keeps track of the actuall rerootings
    rerootings=[]

    #keeps the original tree
    if validRooting(tree, locus_map):
        rerootings+=[tree]
    originaltree=frozenset([get_name(left_subtree(tree)),get_name(right_subtree(tree))])
    alltrees.add(originaltree)

    #call the helper function where actually rerooting happens
    helperRerooting(tree, alltrees, rerootings, locus_map)
    return rerootings  

def helperRerooting(tree, alltrees, rerootings, locus_map):
    """
    The actual function for rerooting
    Reroot on the path to grandchildren of root if possible
    calls addTree to record rerootings and trees generated so far

    Parameters
    --------------------
        tree        -- a 4-tuple tree

    """
    
    if not is_leaf(left_subtree(tree)):
        #if the left tree has two children
        #first we can root on the path to its left children
        newleft=left_subtree(left_subtree(tree))
        newright=(get_name(left_subtree(tree)),right_subtree(left_subtree(tree)),right_subtree(tree))

        #if this tree has not been encountered, add it
        addTree(tree,newleft,newright, alltrees, rerootings, locus_map)
    
        #second we can root on the path to its left
        newleft=(get_name(left_subtree(tree)),left_subtree(left_subtree(tree)),right_subtree(tree))
        newright=right_subtree(left_subtree(tree))

        addTree(tree,newleft,newright, alltrees, rerootings, locus_map)
        

    if not is_leaf(right_subtree(tree)):
        #if the right tree has two children
        #first we can root on the path to its right children
        newright=right_subtree(right_subtree(tree))
        newleft=(get_name(right_subtree(tree)),left_subtree(right_subtree(tree)),left_subtree(tree))

        addTree(tree,newleft,newright, alltrees, rerootings, locus_map)
    
        #second we can root on the path to its left
        newright=(get_name(right_subtree(tree)),right_subtree(right_subtree(tree)),left_subtree(tree))
        newleft=left_subtree(right_subtree(tree))

        addTree(tree,newleft,newright, alltrees, rerootings, locus_map)

def addTree(tree, newleft, newright, alltrees, rerootings, locus_map):
    """
    helper function for the helperRerooting function
    This function add a new rerooting to the rerootings and alltrees if not already found
    and calls the helperRerooting on the new rerooting

    Parameters
    --------------------
        tree        -- a 4-tuple tree that is the orginal tree
        newleft     -- a new left tree generated in helperRerooting
        newright    -- a new right tree generated in helperRerooting

    """
    newtree=frozenset([get_name(newright),get_name(newleft)])
    if newtree not in alltrees:
        newrooting=(get_name(tree),newleft,newright)
        if validRooting(newrooting, locus_map):
            rerootings+=[newrooting]
        alltrees.add(newtree)
        helperRerooting(newrooting, alltrees, rerootings, locus_map)

def parseTreeForDP(tree, parasite):
    """
    Input: four tuple tree as defined by the trees.py (xenoGI)
    Output: tree format represented by a dictionary as required by the 
            DPLOT_DP.py. Note, if an edge directs to a tip, than the children 
            edges will be None
    """
    def parseHelper(tree,treeDict):
        """
        recursive helper. Note the edge to the root of 
        the tree is already in treeDict
        """
        if tree[1]==(): pass  #can't be at a leaf
        lt=tree[1]
        rt=tree[2]
        if lt[1]==(): #if left child is a leaf, add the edge to left child and stop
            key=lt[0]
            value=(tree[0],lt[0],None,None)
            treeDict[key]=value
        else: #if left child is not a leaf, add the edge to left child and recur on lt
            key=lt[0]
            #the left child of left tree
            lt_lt=lt[1]
            rt_lt=lt[2]
            value=(tree[0],lt[0],lt_lt[0],rt_lt[0])
            treeDict[key]=value
            treeDict=parseHelper(lt,treeDict)
        if rt[1]==(): #if right child is a leaf, add the edge to right child and stop
            key=rt[0]
            value=(tree[0],rt[0],None,None)
            treeDict[key]=value
        else: #if right child is not a leaf, add the edge to right child and recur on lt
            key=rt[0]
            lt_rt=rt[1]
            rt_rt=rt[2]
            value=(tree[0],rt[0],lt_rt[0],rt_rt[0])
            treeDict[key]=value
            treeDict=parseHelper(rt,treeDict)

        return treeDict



    if tree[1]==(): #this tree cannot be only a root
        return None
    parsedTree={}
    
    
    lt=tree[1]
    rt=tree[2]
    if parasite:
    # (start vertex, end vertex, left child edge name, right child edge name)
        value=("p_root",tree[0], lt[0], rt[0])
        key="pTop"
    else:
        value=("h_root",tree[0], lt[0], rt[0])
        key="hTop"
    parsedTree[key]=value
    
    return parseHelper(tree, parsedTree)
    
        
        
    
