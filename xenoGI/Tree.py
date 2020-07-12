from Bio import Phylo

## Globals

PARENTOFROOT = "" # name for the parent of the root node in rooted trees

## Classes

class Tree:
    def __init__(self, nodeConnectD=None):
        '''Base class to be inherited by Rtree and Utree.'''
        self.nodeConnectD = nodeConnectD
        self.leafNodeT = None
        self.internalNodeT = None
        self.preOrderT = None

        if self.nodeConnectD != None:
            self.__updateSecondaryAttributes__()

    def preorder(self):
        return self.preOrderT

    def leaves(self):
        return self.leafNodeT

    def internals(self):
        return self.internalNodeT

    def leafCount(self):
        return len(self.leafNodeT)

    def internalNodeCount(self):
        return len(self.internalNodeT)
    
    def nodeCount(self):
        return len(self.nodeConnectD)

    def multifurcatingNodes(self):
        '''Return a list of all nodes where there's a multifurcation (more
than 3 branches).'''
        multifurcL = []
        for node,connecT in self.nodeConnectD.items():
            if len(connecT) > 3:
                multifurcL.append(node)
        return multifurcL

    def isLeaf(self,node):
        '''Return boolean if node is a leaf'''
        if len(self.nodeConnectD[node]) == 1:
            return True
        else: return False

    def children(self,node):
        '''Return tuple of children of node.'''
        connecT = self.nodeConnectD[node]
        return connecT[1:]
        
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

    def __traversePreOrderNodeConnectD__(self,D,node,parentNode):
        '''Traverse nodeConnectD starting at node. Do not recurse along
parentNode.
        '''
        connecT=D[node]
        if len(connecT)==1:
            return [node]
        else:
            outL = [node]
            for child in connecT:
                if child != parentNode:
                    tempL = self.__traversePreOrderNodeConnectD__(D,child,node)
                    outL = outL + tempL
            return outL

    def __traverseForNewickStr__(self,node,parentNode):
        connecT=self.nodeConnectD[node]
        if len(connecT)==1: # tip
            return node
        else:
            outL = []
            for child in connecT:
                if child != parentNode:
                    newickStr = self.__traverseForNewickStr__(child,node)
                    outL = outL + [newickStr]
            return "("+ ",".join(outL)+")"+node

class Rtree(Tree):
    def __init__(self, nodeConnectD=None,rootNode=None):
        '''Initialize a rooted tree object.'''
        super().__init__(nodeConnectD)
        self.rootNode = rootNode

        if self.nodeConnectD != None:
            self.preOrderT = self.__traversePreOrder__(self.rootNode)
        
    def fromString(self,treeStr):
        '''Populate attributes by parsing the string treeStr (which has likely
been read from a file). This method allows multifurcating trees.'''

        L = treeStr.split(",")
        self.rootNode = L[0]
        self.nodeConnectD = {}
        for entryStr in L[1:]:
            entryL = entryStr.split(" ")
            node=entryL[0]
            connecT=tuple(entryL[1:])
            self.nodeConnectD[node] = connecT
        self.__updateSecondaryAttributes__()
            
    def fromNewickFileLoadSpeciesTree(self,treeFN):
        '''Populate attributes based on newick file in treeFN. This method
assumes we are working with a species tree. It must be rooted, bifurcating
and have named internal nodes.
        '''
        bpTree = Phylo.read(treeFN, 'newick', rooted=True)
        self.__checkSpeciesTree__(bpTree)
        self.nodeConnectD = self.__bioPhyloToNodeConnectD__(bpTree)
        self.rootNode = bpTree.root.name
        self.__updateSecondaryAttributes__()
        self.preOrderT = self.__traversePreOrder__(self.rootNode)

    def toNewickStr(self):
        '''Output a newick string.'''
        return self.__traverseForNewickStr__(self.rootNode,PARENTOFROOT)
        
    def subtree(self,node):
        '''Return a new Rtree object with the subtree rooted at node.'''
        
        def traverse(D,subD,node):
            '''Get nodeConnectD for subtree'''
            connecT=D[node]
            subD[node] = connecT
            if len(connecT)==1:
                return
            else:
                for child in connecT[1:]: # 1st entry is parent
                    traverse(D,subD,child)
                return

        subD = {}
        traverse(self.nodeConnectD,subD,node)

        # put special value in parent of root position
        oldConnecT = subD[node]
        newConnecT = (PARENTOFROOT,)+oldConnecT[1:]
        subD[node] = newConnecT
        
        #return Rtree(subD,node)
        return subD
        
    def createSubtreeD(self):
        '''Get all subtrees and put them in a dict keyed by node name.'''
        subtreeD = {}
        for node in self.preorder():
            subtree = self.subtree(node)
            subtreeD[node] = subtree
        return subtreeD
            
    def ancestors(self,node):
        '''Return a tuple of nodes ancestral to node.'''

        def traverse(D,node):
            connecT=D[node]
            parent = connecT[0]
            if parent == PARENTOFROOT: # at root
                return []
            else:
                return [parent]+traverse(D,parent)

        return tuple(traverse(self.nodeConnectD,node))

    def getParent(self,node):
        connecT = self.nodeConnectD[node]
        return connecT[0] # in rooted tree, parent is first element
        
    def getNearestNeighborL(self,leaf):
        '''Given a leaf, return a list containing the other leaf or
    leaves which are most closely related to leaf.'''
        parent = self.getParent(leaf) # assume leaf is in tree
        subRtree = self.subtree(parent)
        leafL = list(subRtree.leaves())
        leafL.remove(leaf)
        return leafL

    def findMrcaPair(self,node1,node2):
        '''Return the most recent common ancestor of node1 and node2.'''

        anc1L = (node1,) + self.ancestors(node1)
        anc2L = (node2,) + self.ancestors(node2)

        for node in anc1L:
            if node in anc2L:
                return node
        return None

    def findMrca(self,nodeL):
        '''Get mrca of a list of nodes'''
        mrca = nodeL[0]
        for node in nodeL[1:]:
            mrca = self.findMrcaPair(mrca,node)
        return mrca
            
    def createDtlorD(self,spTreeB):
        '''Return a dict in the format expected by the dtlor reconciliation
code. Key is a node, value is (parent, node, child1, child2). For now
this only works for bifurcating trees. spTreeB is a boolean that is
True if this is a species tree, and False if it is for a gene
tree.
        '''
        if self.multifurcatingNodes() != []:
            raise ValueError("This tree is not bifurcating.")
        
        dtlorD = {}
        for node,connecT in self.nodeConnectD.items():
            if len(connecT) == 1:
                # leaf
                dtlorD[node] = (parent, node, None, None)
            else:
                parent,child1,child2 = connecT                
                if node == self.rootNode:
                    if spTreeB:
                        dtlorD[node] = ("h_root", node, child1, child2)
                    else:
                        dtlorD[node] = ("p_root", node, child1, child2)
                else:
                    dtlorD[node] = (parent, node, child1, child2)
                    
        return dtlorD
        
    def fileStr(self):
        '''Return a string representation of an Rtree object. This consists of
the root node, followed by a series of entrys from nodeConnectD. These
elements are separated by commas (and the entries in nodeConnectD are
internally separated by spaces.
        '''
        outL=[self.rootNode] # first thing is root node.
        for node,connecT in self.nodeConnectD.items():
            # elements in one entry of nodeConnectD separated by spaces
            entry = node + " " + " ".join(connecT)
            outL.append(entry)
        return ",".join(outL)
            
    ## methods not for end users
    
    def __checkSpeciesTree__(self,bpTree):
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
        
    def __bioPhyloToNodeConnectD__(self,bpTree):
        '''Convert a biopython tree object to a node connection dict with keys
    that are nodes, and values that are (parent, child1, child2...).'''
        nodeConnectD = {}
        # parent of root is PARENTOFROOT
        return self.__bioPhyloCladeToNodeConnectD__(bpTree.root,nodeConnectD,PARENTOFROOT)

    def __bioPhyloCladeToNodeConnectD__(self,clade,nodeConnectD,parent):
        '''Recursive helper to convert a biopython clade object.
        '''
        if clade.is_terminal():
            nodeConnectD[clade.name] = (parent,)
            return nodeConnectD
        else:
            nodeConnectD[clade.name] = (parent,clade[0].name,clade[1].name)
            for iterClade in clade:
                nodeConnectD = self.__bioPhyloCladeToNodeConnectD__(iterClade,nodeConnectD,clade.name)
            return nodeConnectD

    def __traversePreOrder__(self,node):
        '''Traverse in preorder starting at rootNode'''
        return tuple(self.__traversePreOrderNodeConnectD__(self.nodeConnectD,self.rootNode,PARENTOFROOT))
    def __repr__(self):
        return "Rtree: "+self.toNewickStr()


class Utree(Tree):
    def __init__(self, nodeConnectD=None):
        '''Initialize an unrooted tree object.'''
        super().__init__(nodeConnectD)
        
        if self.nodeConnectD != None:
            # arbitraryNode provides a starting point for some methods
            self.arbitraryNode = self.internals()[0]
            self.preOrderT = self.__traversePreOrder__(self.arbitraryNode)
            self.branchPairT = self.__createBranchPairT__()
            
    def fromNewickFile(self,treeFN):
        '''Populate attributes based on newick file in treeFN. This method
assumes we are working with an unrooted gene tree that does not have
named interal nodes (we will create those).

        '''
        bpTree = Phylo.read(treeFN, 'newick', rooted=False)
        iNodeNum,nodeConnectD = self.__bioPhyloToNodeConnectD__(bpTree.clade,PARENTOFROOT,{},0)
        self.nodeConnectD = nodeConnectD
        self.__updateSecondaryAttributes__()
        self.arbitraryNode = self.internals()[0]
        self.preOrderT = self.__traversePreOrder__(self.arbitraryNode)
        self.branchPairT = self.__createBranchPairT__()

    def toNewickStr(self):
        '''Output a newick string.'''
        return self.__traverseForNewickStr__(self.arbitraryNode,PARENTOFROOT)

    def root(self,branchPair):
        '''Root using the tuple branchPair and return an Rtree object.'''

        def traverse(D,newD,node,parentNode):
            '''Get nodeConnectD for the part of the tree defined by node, and in
the oposite direction from parentNode.
            '''
            oldConnecT=D[node]

            # make sure parent node is first
            assert(parentNode in oldConnecT)
            connecL = [parentNode]
            for tempNode in oldConnecT:
                if tempNode != parentNode:
                    connecL.append(tempNode)
            connecT = tuple(connecL)
            
            newD[node] = oldConnecT # store
            
            if len(connecT)==1:
                return
            else:
                for child in connecT:
                    if child != parentNode:
                        traverse(D,newD,child,node)
                return

        def updateChildOfRoot(newD,nodeToWorkOn,nodeToReplace):
            '''Adjust child of root to say root is parent.'''
            newConnecL = ["root"]
            for node in newD[nodeToWorkOn]:
                if node != nodeToReplace:
                    newConnecL.append(node)
            newD[nodeToWorkOn] = tuple(newConnecL)
            return newD
            
        # make new nodeConnectD
        newD = {}
        traverse(self.nodeConnectD,newD,branchPair[0],branchPair[1])
        traverse(self.nodeConnectD,newD,branchPair[1],branchPair[0])

        # must add the root
        rootNode = "root"
        newD[rootNode] = (PARENTOFROOT,branchPair[0],branchPair[1])

        # adjust children of root to say root is parent
        newD = updateChildOfRoot(newD,branchPair[0],branchPair[1])
        newD = updateChildOfRoot(newD,branchPair[1],branchPair[0])
        
        return Rtree(newD,rootNode)

    def iterAllRootedTrees(self):
        '''Iterator yielding all possible rooted trees from this unrooted tree.'''
        for branchPair in self.branchPairT:
            yield self.root(branchPair)
        
    def __bioPhyloToNodeConnectD__(self,bpClade,parentNodeStr,nodeConnectD,iNodeNum):
        '''Convert a biopython clade object to a node connection dict with
keys that are nodes, and values that are the nodes connected to. For
internals (node1, node2, node3...). We assume bpClade does not have
named internal nodes, and we name them here: g0, g1 etc.

        '''
        if parentNodeStr == PARENTOFROOT:
            connecL = [] # for the outer case, where we're not coming from anywhere
        else:
            connecL = [parentNodeStr]

        if bpClade.is_terminal():
            thisNodeStr = bpClade.name
            nodeConnectD[thisNodeStr] = tuple(connecL)
        else:
            thisNodeStr = "g"+str(iNodeNum)
            iNodeNum += 1
            for childClade in bpClade:

                # collect connection
                if childClade.is_terminal():
                    connecL.append(childClade.name)
                else:
                    # the connecting node will be the next one (iNodeNum)
                    connecL.append("g"+str(iNodeNum))

                # recurse
                iNodeNum,nodeConnectD = self.__bioPhyloToNodeConnectD__(childClade,thisNodeStr,nodeConnectD,iNodeNum)

            # add in thie connection for this node
            nodeConnectD[thisNodeStr] = tuple(connecL)

        return iNodeNum,nodeConnectD

    def __createBranchPairT__(self):
        '''Make the branchPairT attribute to represent the branches. This is a
tuple of tuples, where the subtuples represent edges and consist of
the two nodes on either end.
        '''
        S = set()
        for node,connecT in self.nodeConnectD.items():
            for otherNode in connecT:
                # put them in preorder
                edgeL=[(node,self.preOrderT.index(node)),(otherNode,self.preOrderT.index(otherNode))]                
                edgeL.sort(key=lambda x: x[1])
                edgeT = tuple(x[0] for x in edgeL)
                S.add(edgeT)
        return tuple(sorted(S))
        
    def __traversePreOrder__(self,node):
        '''Traverse in preorder starting at arbitraryNode'''
        return tuple(self.__traversePreOrderNodeConnectD__(self.nodeConnectD,self.arbitraryNode,PARENTOFROOT))
    def __repr__(self):
        return "Utree: "+self.toNewickStr()
