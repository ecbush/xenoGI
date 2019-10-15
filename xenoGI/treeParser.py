#file for parsing the tree into a format tha works for the DTLOR algorithm
# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name (will use end vertex to represent) and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop". 
# Similarly, the dummy edge leading to the root of the host tree must be named
#"hTop"


# tree=("root",("left",("HI",(),(),2),("HA",(),(),1),1),('right',(),(),2),3)
def parseTreeForDP(tree, parasite=True):
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
    
    key=tree[0]
    lt=tree[1]
    rt=tree[2]
    if parasite:
    # (start vertex, end vertex, left child edge name, right child edge name)
        value=("pTop",tree[0], lt[0], rt[0])
    else:
        value=("hTop",tree[0], lt[0], rt[0])
    parsedTree[key]=value
    
    return parseHelper(tree, parsedTree)
    
        
        
    
