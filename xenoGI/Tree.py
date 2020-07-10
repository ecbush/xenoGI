
from Bio import Phylo


class Rtree:
    def __init__(self, nodeConnectD=None,rootNode=None):
        '''Initialize a rooted tree object.'''
        self.nodeConnectD = nodeConnectD
        self.rootNode = rootNode
        self.leafNodeT = None
        self.internalNodeT = None
        self.preOrderT = None
        
        if self.nodeConnectD != None:
            self.__updateSecondaryAttributes__()
        
    def __updateSecondaryAttributes__(self):
        '''Create leafNodeT and internalNodeT.'''
        leafNodeL = []
        internalNodeL = []
        for node,connecT in self.nodeConnectD.items():
            if len(connecT) == 1:
                leafNodeL.append(node)
            else:
                internalNodeL.append(node)
        self.leafNodeT = tuple(leafNodeL)
        self.internalNodeT = tuple(internalNodeL)
        self.preOrderT = self.__traversePreOrder__()
        
    def readNewickSpeciesTree(self,treeFN):
        '''Fill in our attributes based on newick file in treeFN. This tree
must be rooted, bifurcating and have named internal nodes.'''

        # read biopython
        bpTree = Phylo.read(treeFN, 'newick', rooted=True)
        checkSpeciesTree(bpTree)
        self.nodeConnectD = bioPhyloToNodeConnectD(bpTree)
        self.rootNode = bpTree.root.name
        self.__updateSecondaryAttributes__()

    def __traversePreOrder__(self):
        # support func
        def traverse(D,node):
            connecT=D[node]
            if len(connecT)==1:
                return [node]
            else:
                outL = [node]
                for child in connecT[1:]:
                    tempL = traverse(D,child)
                    outL = outL + tempL
                return outL
        # end support func
        return tuple(traverse(self.nodeConnectD,self.rootNode))

    def preorder(self):
        return self.preOrderT

    def leaves(self):
        return self.leafNodeT

    def internals(self):
        return self.internalNodeT

    def leafCount(self):
        return len(self.leafNodeT)

    def nodeCount(self):
        return len(self.nodeConnectD)

    def subtree(self,node):
        '''Return a new Tree object with the subtree rooted at node.'''
        # get nodeConnectD for subtree
        def traverse(D,node,subD):
            connecT=D[node]
            subD[node] = connecT
            if len(connecT)==1:
                return
            else:
                for child in connecT[1:]:
                    traverse(D,child,subD)
                return
        # end support

        subD = {}
        traverse(self.nodeConnectD,node,subD)

        return Rtree(subD,node)
    
## Biopython funcs
def bioPhyloToNodeConnectD(bpTree):
    '''Convert a biopython tree object to a node connection dict with keys
that are nodes, and values that are (parent, child1, child2...).'''
    nodeConnectD = {}
    return bioPhyloCladeToNodeConnectD(bpTree.root,nodeConnectD,None)

def bioPhyloCladeToNodeConnectD(clade,nodeConnectD,parent):
    '''Recursive helper to convert a biopython clade object. We can assume
tree is bifurcating

    '''
    if clade.is_terminal():
        nodeConnectD[clade.name] = (parent,)
        return nodeConnectD
    else:
        nodeConnectD[clade.name] = (parent,clade[0].name,clade[1].name)
        nodeConnectD = bioPhyloCladeToNodeConnectD(clade[0],nodeConnectD,clade.name)
        nodeConnectD = bioPhyloCladeToNodeConnectD(clade[1],nodeConnectD,clade.name)
        return nodeConnectD

def checkSpeciesTree(bpTree):
    '''Check that a biopython tree is rooted, bifurcating, and has named
internal nodes. Throw error if not. Returns None.

    '''

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


## main
#bpTree = Phylo.read("testE.tre", 'newick', rooted=True)
#tupleTree=bioPhyloToTupleTree(bpTree)

rtree = Rtree()
rtree.readNewickSpeciesTree("testE.tre")
