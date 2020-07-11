from Bio import Phylo

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

    def nodeCount(self):
        return len(self.nodeConnectD)

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

    def __traversePreOrder__(self,node):
        
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
        
        return tuple(traverse(self.nodeConnectD,node))

    def __traverseForNewickStr__(self,node):
        connecT=self.nodeConnectD[node]
        if len(connecT)==1: # tip
            return node
        else:
            outL = []
            for child in connecT[1:]:
                newickStr = self.__traverseForNewickStr__(child)
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

    def toNewickStr(self):
        '''Output a newick string.'''
        return self.__traverseForNewickStr__(self.rootNode)
        
    def subtree(self,node):
        '''Return a new Tree object with the subtree rooted at node.'''
        
        def traverse(D,node,subD):
            '''Get nodeConnectD for subtree'''
            connecT=D[node]
            subD[node] = connecT
            if len(connecT)==1:
                return
            else:
                for child in connecT[1:]:
                    traverse(D,child,subD)
                return

        subD = {}
        traverse(self.nodeConnectD,node,subD)

        return Rtree(subD,node)

    def ancestors(self,node):
        '''Return a tuple of nodes ancestral to node.'''

        def traverse(D,node):
            connecT=D[node]
            parent = connecT[0]
            if parent == "": # at root
                return []
            else:
                return [parent]+traverse(D,parent)

        return tuple(traverse(self.nodeConnectD,node))

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
            
    ## methods that not for end users
    
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
        # parent of root is ""
        return self.__bioPhyloCladeToNodeConnectD__(bpTree.root,nodeConnectD,"")

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

    def __repr__(self):
        return "Rtree: "+self.toNewickStr()


class Utree(Tree):
    def __init__(self, nodeConnectD=None):
        '''Initialize an unrooted tree object.'''
        super().__init__(nodeConnectD)
        
        if self.nodeConnectD != None:
            self.arbitraryNode = self.internals()[0]
            self.preOrderT = self.__traversePreOrder__(self.arbitraryNode)

    def fromNewickFile(self,treeFN):
        '''Populate attributes based on newick file in treeFN. This method
assumes we are working with an unrooted gene tree that does not have
named interal nodes (we will create those).

        '''
        bpTree = Phylo.read(treeFN, 'newick', rooted=False)
        self.nodeConnectD = {}
        iNodeNum = 0
        iNodeNum = self.__bioPhyloToNodeConnectD__(bpTree.clade,iNodeNum,"g"+str(iNodeNum))
        self.__updateSecondaryAttributes__()
        self.arbitraryNode = self.internals()[0]
        self.preOrderT = self.__traversePreOrder__(self.arbitraryNode)
 
    def toNewickStr(self):
        '''Output a newick string.'''
        return self.__traverseForNewickStr__(self.arbitraryNode)

    def __bioPhyloToNodeConnectD__(self,bpClade,iNodeNum,parentNodeStr):
        '''Convert a biopython clade object to a node connection dict with keys
    that are nodes, and values that are the nodes connected to. For internals (node1, node2, node3...).'''

        iNodeNum += 1
        thisNodeStr = "g"+str(iNodeNum)
        connecL = [parentNodeStr]

        if bpClade.is_terminal():
            self.nodeConnectD[thisNodeStr] = tuple(connecL)
        else:
            for childClade in bpClade:

                # collect connection
                if childClade.is_terminal():
                    connecL.append(childClade.name)
                else:
                    # the connecting node will be the next one. but don't increment yet
                    connecL.append("g"+str(iNodeNum+1))

                # recurse
                iNodeNum = self.__bioPhyloToNodeConnectD__(childClade,iNodeNum,thisNodeStr)

            # add in thie connection for this node
            self.nodeConnectD[thisNodeStr] = tuple(connecL)
            
        return iNodeNum
