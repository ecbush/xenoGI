# functions for loading and manipulating phylogenetic trees
from Bio import Phylo

def bioPhyloToTupleTree(bpTree):
    '''Convert a biopython tree object to 4 tuple tree.'''
    return bioPhyloCladeToTupleTree(bpTree.root)

def bioPhyloCladeToTupleTree(clade):
    '''Convert a biopython clade object to 4 tuple tree.'''
    nm = clade.name
    br = clade.branch_length
    if clade.is_terminal():
        return (nm,(),(),br)
    else:
        lt = bioPhyloCladeToTupleTree(clade[0])
        rt = bioPhyloCladeToTupleTree(clade[1])
        return (nm,lt,rt,br)
    
def nodeCount(tree):
    '''How many nodes in tree?'''
    if tree[1]==():
        return 1
    else:
        return 1+ nodeCount(tree[1]) + nodeCount(tree[2])

def nodeList(tree):
    '''Return list of nodes in tree.'''
    if tree[1]==():
        return [tree[0]]
    else:
        return [tree[0]] + nodeList(tree[1]) + nodeList(tree[2])
    
def leafCount(tree):
    '''How many leaves in tree?'''
    if tree[1]==():
        return 1
    else:
        return leafCount(tree[1]) + leafCount(tree[2])

def leafList(tree):
    '''Return list of leaves in tree.'''
    if tree[1]==():
        return [tree[0]]
    else:
        return leafList(tree[1]) + leafList(tree[2])

def iNodeList(tree):
    '''Return list of internal nodes in tree.'''
    if tree[1]==():
        return []
    else:
        return [tree[0]] + iNodeList(tree[1]) + iNodeList(tree[2])
    
def subtree(tree,node):
    '''Return the subtree with node at its root. Assume node is in tree.'''
    if tree[0]==node:
        return tree
    elif tree[1]==():
        return None
    else:
        l=subtree(tree[1],node)
        r=subtree(tree[2],node)
        if l==None:
            return r
        else:
            return l
    
def isRootNode(tree,mrcaNum):
    '''Is mrcaNum the root node?'''
    # root node is tree[0] in our tuple trees
    return mrcaNum == tree[0]
    
def strTree2numTree(tree,counter):
    '''Given a tuple tree with nodes specified by strings, convert to
numbered nodes. Return tree with numbered nodes.'''
    if tree[1]==():
        return (counter,(),(),tree[3]),counter+1
    else:
        leftNumTree,counter=strTree2numTree(tree[1],counter)
        rightNumTree,counter=strTree2numTree(tree[2],counter)
        numTree=(counter,leftNumTree,rightNumTree,tree[3])
        return numTree,counter+1

def checkTree(bpTree):
    '''Check that a biopython tree is rooted and has named internal
nodes. Throw error if not. Returns None.'''

    # check its bifurcating. This actually ignors the root, but checks
    # for non-bifurcating nodes below that.
    if not bpTree.is_bifurcating():
        raise ValueError("This tree is not bifurcating.")
    
    # Now check to make sure there are only two clades at the root of
    # this tree. Otherwise, its going to mess up bioPhyloToTupleTree
    if len(bpTree.clade) > 2:
        raise ValueError("readTree requires a rooted tree with exactly two clades at the root. This input tree has more than two. This typically means that it is intended to be read as an unrooted tree.")

    # Make sure the internal nodes are labelled.
    predictedNumIntNodes = bpTree.count_terminals() - 1
    internalL = bpTree.get_nonterminals()
    if len(internalL) != predictedNumIntNodes or any((n.name==None for n in internalL)):
        raise ValueError("All the internal nodes of the input tree must have names.")
    
def readTree(filename):
    '''Read Newick tree from file and convert it to tuple format.'''

    bpTree = Phylo.read(filename, 'newick', rooted=True)

    checkTree(bpTree)
    
    stringTree=bioPhyloToTupleTree(bpTree)
    counter=0
    numTree,counter=strTree2numTree(stringTree,counter)

    strainNamesO = strainNames(stringTree,numTree)
    
    return numTree,strainNamesO

def createSubtreeL(tree):
    '''Return a list containing all subtrees.'''
    if tree[1]==():
        return [tree]
    else:
        return [tree]+createSubtreeL(tree[1]) + createSubtreeL(tree[2])

def getNearestNeighborL(leaf,tree):
    '''Given a leaf and a tree, return a list containing the other leaf or
leaves which are most closely related to leaf.'''

    parent = getParent(leaf,tree) # assume leaf is in tree
    leafL = leafList(subtree(tree,parent))
    leafL.remove(leaf)
    return leafL
    
def getParent(leaf,tree):
    '''Return parent node of leaf.'''
    if tree[1] == ():
        # We've come too far, it wasn't there
        return None
    elif tree[1][0] == leaf or tree[2][0] == leaf:
        return tree[0]
    else:
        left = getParent(leaf,tree[1])
        right = getParent(leaf,tree[2])

        if left != None:
            return left
        else: return right
        
def tupleTree2Newick(tree):
    '''Convert a four tuple based tree (root,left,right,branchLen) into a
newick formated string.'''
    if tree[1]==():
        return str(tree[0])+":"+str(tree[3])
    else:
        leftStr=tupleTree2Newick(tree[1])
        rightStr=tupleTree2Newick(tree[2])
        return "("+leftStr+","+rightStr+"):"+str(tree[3])

def writeTree(tree,fileName):
    '''Write tree to fileName (in newick format).'''
    f=open(fileName,"w")
    f.write("("+tupleTree2Newick(tree)+");")
    f.close()

class strainNames:

    def __init__(self,stringTree,numTree):
        '''Create an object for converting between the numerical and string
versions of strain names.'''

        # make dictionaries for converting between number and string strain names
        self.nameToNumD={}
        self.makeTreeD(stringTree,numTree,self.nameToNumD)
        self.numToNameD={}
        self.makeTreeD(numTree,stringTree,self.numToNameD)

    def makeTreeD(self,tree1,tree2,treeD):
        '''Make a dictionary to convert from node names in tree1 to node names
    in the identically shaped tree2.'''
        treeD[tree1[0]]=tree2[0]
        if tree1[1]==(): return
        else:
            self.makeTreeD(tree1[1],tree2[1],treeD)
            self.makeTreeD(tree1[2],tree2[2],treeD)
            return

    def nameToNum(self,strainStr):
        return self.nameToNumD[strainStr]

    def numToName(self,strainNum):
        return self.numToNameD[strainNum]

    def isStrainPresentByName(self,strainName):
        return strainName in self.nameToNumD

    def isStrainPresentByNum(self,strainNum):
        return strainNum in self.numToNameD

    def iterStrainNums(self):
        for strainNum in self.numToNameD.keys():
            yield strainNum

    def iterStrainNames(self):
        for strainName in self.nameToNumD.keys():
            yield strainName
