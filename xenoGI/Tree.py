from Bio import Phylo
from collections import OrderedDict
from xenoGI import trees

## Globals

ROOT_PARENT_NAME = "" # name for the parent of the root node in rooted trees

## Classes

class Tree:
    def __init__(self, nodeConnectD=None):
        '''Base class to be inherited by Rtree and Utree. Key internal data
structure is nodeConnectD, which is keyed by node and has values
(parent, child1, child2...). All nodes names are represented as
strings. In the case of a gene tree, the tip may be a gene number, but
we still hold it as a string.

        '''
        self.nodeConnectD = nodeConnectD
        self.leafNodeT = None
        self.internalNodeT = None
        self.preOrderT = None
        self.branchLenD = None
        
        if self.nodeConnectD != None:
            self.__updateSecondaryAttributes__()

    def fileStr(self):
        '''Return a string representation of a Tree object. Outer separator is "|".

        rootNode/arbitraryNode | node connection info | branch len info

        The first position gives the root node in a rooted tree, or
        else the arbitraryNode in an unrooted tree.

        Node connection info is organized as follows. Each entry in
        nodeConnectD gets internally separated by spaces. Different
        entries are separated by commas.

        Branch connection info is organized similarly. Each entry in
        branchLenD gets internally separated by spaces. Different
        entries are separated by commas. If branch len info
        is absent, will be None.
        '''
        outL = []

        # root/arbitary node
        if hasattr(self,"rootNode"):
            # rooted tree
            outL.append(str(self.rootNode))
        else:
            # unrooted
            outL.append(str(self.arbitraryNode))

        # node connection
        nodeConnecL=[] # first thing is root node.
        for node,connecT in self.nodeConnectD.items():
            # elements in one entry of nodeConnectD separated by spaces
            entry = node + " " + " ".join(connecT)
            nodeConnecL.append(entry)
        nodeConnecStr = ",".join(nodeConnecL)
        outL.append(nodeConnecStr)

        # branches
        if self.branchLenD == None:
            outL.append("None")
        else:
            branchLenL = []
            for branchPair,brLen in self.branchLenD.items():
                entry = branchPair[0] + " " + branchPair[1] + " " + str(brLen)
                branchLenL.append(entry)
            branchLenStr = ",".join(branchLenL)
            outL.append(branchLenStr)
            
        return "|".join(outL)

    def populateAttributes(self,nodeConnectD,rootArbNode,branchLenD=None):
        '''Populate the attributes of an empty tree object.'''

        if self.nodeConnectD != None or self.branchLenD != None:
            raise ValueError("Attempting to populate attributes of a non-empty tree.")    

        self.nodeConnectD = nodeConnectD
        self.__updateSecondaryAttributes__()
        self.preOrderT = self.__traversePreOrder__(rootArbNode)
        
        # root/arbitraryNode
        if isinstance(self, Rtree):
            self.rootNode = rootArbNode
        else:
            self.arbitraryNode = rootArbNode
            
        # branches
        if branchLenD == None:
            self.branchLenD = branchLenD
        else:
            self.branchPairT = self.__createBranchPairT__()
            self.branchLenD = {}
            for branchPair in self.branchPairT:
                if branchPair in branchLenD:
                    self.branchLenD[branchPair] = branchLenD[branchPair]
                else:
                    # flip order
                    self.branchLenD[branchPair] = branchLenD[(branchPair[1],branchPair[0])]
        
    def fromString(self,treeStr):
        '''Populate attributes by parsing the string treeStr (which has likely
been read from a file, and is our native file format for trees). This
method allows multifurcating trees. See fileStr method for description
of format.
        '''
        rootArbNode,nodeConnecStr,branchLenStr = treeStr.split("|")

        # nodeConnectD
        L = nodeConnecStr.split(",")
        nodeConnectD = {}
        for entryStr in L:
            entryL = entryStr.split(" ")
            node=entryL[0]
            connecT=tuple(entryL[1:])
            nodeConnectD[node] = connecT

        # branchLenD
        if branchLenStr == "None":
            branchLenD = None
        else:
            L = branchLenStr.split(",")
            branchLenD = {}
            for entryStr in L:
                br0,br1,brLen = entryStr.split(" ")
                branchLenD[(br0,br1)] = float(brLen)

        # put in these attributes
        self.populateAttributes(nodeConnectD,rootArbNode,branchLenD)
                
    def isLeaf(self,node):
        '''Return boolean if node is a leaf'''
        if len(self.nodeConnectD[node]) < 3:
            # internal node must have at least 3 connections
            return True
        else: return False

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

    def getParent(self,node):
        '''For rooted trees, give parent. For unrooted trees, this simply
gives a node leading to this node. We have this method also for
unrooted trees since it can be useful for finding what leads to a
tip.
        '''
        connecT = self.nodeConnectD[node]
        return connecT[0] # parent is first element

    def getBranchLen(self,branchPair):
        '''Return the branch length for the branch defined by branchPair.'''

        if self.branchLenD == None:
            return None
        else:
            return self.branchLenD[branchPair]
        
    def binarize(self,gtLocusMapD=None):
        '''Arbitarily convert a tree with multifurcating nodes into one with
binary nodes. Newly introduced nodes have a "b" in their naming to
mark them off. If tree has branch lengths, then newly introduced
branches are given 0 length. gtLocusMapD gives the syntenic regions of
genes. If this is not None, then we use it to help order the way we
binarize multifurcating nodes.
        '''
        
        # funcs
        def getNewNodes(newNodeCtr,parent,node,childL,updateParentD):
            '''Goes through a list of children from a multifurcating node, and
creates a series of clades by having one group branch off at each
step. First groups listed will be first to branch off. Returns a list
of new nodes, and also an updateParentD which can later be used to
make sure the parents are correct.
            '''
            if len(childL) == 2:
                connecT = (parent,childL[0],childL[1])
                return newNodeCtr,updateParentD,[(node,connecT)],
            else:
                newNode = "b"+str(newNodeCtr)
                newNodeCtr+=1
                connecT = (parent,childL[0],newNode)
                # put rest of children in updateParentD
                for child in childL[1:]:
                    updateParentD[child] = newNode # for now they descent from newNode
                # recurse
                newNodeCtr,updateParentD,L=getNewNodes(newNodeCtr,node,newNode,childL[1:],updateParentD)
                return newNodeCtr,updateParentD,[(node,connecT)]+L

        def orderChildrenBySynteny(childT,gtLocusMapD):
            '''Given a tuple of children from a multifurcating node, and a dict of
 syntenic regions for these, reorder the children in a way that makes
 use of synteny. When we construct the binary tree, the first listed
 children in our output will branch off first. If there are branches
 in childT where synteny is the same in all genes, and where this
 syteny is not shared with other child branches, then we'd like to put
 those first. This tweak is intented to make the binarization more
 likely to reflect what really happened.
            '''
            childToSynRegTupleD = {}
            for child in childT:
                # get the synReg for this child
                if self.isLeaf(child):
                    synRegT = (gtLocusMapD[child],)
                else:
                    # internal node
                    synRegS = set()
                    for chLeaf in self.subtree(child).leaves():
                        synRegS.add(gtLocusMapD[chLeaf])

                    synRegT = tuple(sorted(synRegS))
                    
                childToSynRegTupleD[child] = synRegT

            # create dict keyed by synReg
            synRegToChildD = {}
            for child in childToSynRegTupleD:
                for synReg in childToSynRegTupleD[child]:
                    if synReg in synRegToChildD:
                        synRegToChildD[synReg].append(child)
                    else:
                        synRegToChildD[synReg] = [child]

            # for each child, characterize the max num children its
            # synRegs appear in
            outChildL = []
            for child in childT:

                synRegT = childToSynRegTupleD[child]
                maxAppearances = 0
                for synReg in synRegT:
                    if len(synRegToChildD[synReg]) > maxAppearances:
                        maxAppearances = len(synRegToChildD[synReg])
                
                outChildL.append((maxAppearances,child))
                
            # sort by max appearances, so the ones with less (ie 1) come first
            outChildL.sort(key = lambda x: x[0])

            return tuple((child for _,child in outChildL))
            
        # main part of binarize

        # get node connections
        newNodeConnectD = {}
        updateParentD = {}
        newNodeCtr = 0
        for node,connecT in self.nodeConnectD.items():
            if len(connecT) <= 3:
                # not multifurcating
                newNodeConnectD[node] = connecT
            else:
                # multifurcating        
                parent = connecT[0]
                childT = connecT[1:]

                if gtLocusMapD != None:
                    childT = orderChildrenBySynteny(childT,gtLocusMapD)
                
                newNodeCtr,updateParentD,newNodeL = getNewNodes(newNodeCtr,parent,node,childT,updateParentD)
                for nNode,nConnecT in newNodeL:
                    newNodeConnectD[nNode] = nConnecT

        # fix newNodeConnectD using updateParentD
        for node in updateParentD.keys():
            connecT = newNodeConnectD[node]
            newConnecT = (updateParentD[node],) + connecT[1:] # new parent
            newNodeConnectD[node] = newConnecT
            
        # get branch lens
        if self.branchLenD == None:
            newBranchLenD = None
        else:
            newBranchLenD = {}
            for nNode,nConnecT in newNodeConnectD.items():
                for otherEndOfBranch in nConnecT:
                    if (nNode,otherEndOfBranch) in self.branchLenD:
                        newBranchLenD[(nNode,otherEndOfBranch)] = self.branchLenD[(nNode,otherEndOfBranch)]
                    elif (otherEndOfBranch,nNode) in self.branchLenD:
                        # it could be either order
                        newBranchLenD[(otherEndOfBranch,nNode)] = self.branchLenD[(otherEndOfBranch,nNode)]
                    else:
                        newBranchLenD[(nNode,otherEndOfBranch)] = 0
            
                                
        # create output tree
        if hasattr(self,"rootNode"):      
            # we're a rooted tree
            outTreeO = Rtree()
            rootArbNode = self.rootNode
        else:
            outTreeO = Utree()
            rootArbNode = self.arbitraryNode

        outTreeO.populateAttributes(newNodeConnectD,rootArbNode,newBranchLenD)
        return outTreeO
        
    def __contains__(self,node):
        return node in self.nodeConnectD
        
    def __updateSecondaryAttributes__(self):
        '''Create leafNodeT and internalNodeT.'''
        leafNodeT,internalNodeT = self.__updateSecondaryAttributesHelper__(self.nodeConnectD)
        self.leafNodeT = leafNodeT
        self.internalNodeT = internalNodeT

    def __updateSecondaryAttributesHelper__(self,nodeConnectD):
        '''Actually does the work of getting internals and leaves..'''
        leafNodeL = []
        internalNodeL = []
        for node,connecT in nodeConnectD.items():
            if len(connecT) < 3:
                # it's a tip
                leafNodeL.append(node)
            else:
                internalNodeL.append(node)
        return tuple(leafNodeL),tuple(internalNodeL)
        
    def iterBranches(self):
        '''Iterate over all branches, returning branchPair tuple and branch
length.'''
        if self.branchLenD == None:
            raise ValueError("This unrooted tree object does not have branch lengths defined.")
            
        for branchPair in self.branchPairT:
            branchLen = self.branchLenD[branchPair]
            yield branchPair,branchLen

    def getBranchesByLengthL(self):
        '''Returns a list of tuples (branchPair,branchLen) sorted by branch length. Largest first.'''
        if self.branchLenD == None:
            raise ValueError("This unrooted tree object does not have branch lengths defined.")
            
        L=list(self.branchLenD.items())
        L.sort(key=lambda x: x[1],reverse=True)
        return L
            
    def maxBranchLen(self):
        '''Find the maximum branch length present, and return the
corresponding branchPair tuple and the length.

        '''
        L=list(self.branchLenD.items())
        L.sort(key=lambda x: x[1])
        return L[-1]
        
    def __traversePreOrder__(self,node):
        '''Traverse in preorder starting at node'''
        if self.nodeCount() == 2:
            # two tip tree is a special case            
            return self.leafNodeT
        else:
            return tuple(self.__traversePreOrderNodeConnectD__(self.nodeConnectD,node,ROOT_PARENT_NAME))

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

    def __traverseForNewickStr__(self,node,parentNode,includeBrLength,nodeLabelD):
        '''Produce newick representation of self starting at node. If
includeBrLength is True, include branch lengths. nodeLabelD is a
dictionary containing labels to be added to the nodes. If not to be
used, these should be set to None.
        '''
        connecT=self.nodeConnectD[node]

        # construct node label
        if nodeLabelD == None:
            nodeLab = node
        else:
            nodeLab = node+"["+nodeLabelD[node]+"]"
            
        # construct branch label
        if includeBrLength and parentNode != ROOT_PARENT_NAME:
            brLab = ":"+str(self.branchLenD[(parentNode,node)])
        else:
            brLab = ""
        
        # construct newick    
        if len(connecT)==1: # tip
            return nodeLab+brLab
        else:
            outL = []
            for child in connecT:
                if child != parentNode:
                    newickStr = self.__traverseForNewickStr__(child,node,includeBrLength,nodeLabelD)
                    outL.append(newickStr)
            return "("+ ",".join(outL)+")"+nodeLab+brLab
        
    def __createBranchPairT__(self):
        '''Make the branchPairT attribute to represent the branches. This is a
tuple of tuples, where the subtuples represent edges and consist of
the two nodes on either end.
        '''
        S = set()
        for node,connecT in self.nodeConnectD.items():
            for otherNode in connecT:
                # don't include any edge with ROOT_PARENT_NAME
                if otherNode != ROOT_PARENT_NAME:
                    # put them in preorder
                    edgeT = tuple((nd for nd in self.preorder() if nd in [node,otherNode]))
                    S.add(edgeT)
        return tuple(sorted(S))
        
    def __eq__(self,other):
        '''See if two trees are equivalent.'''
        if hasattr(self,"rootNode") and  hasattr(other,"rootNode"):
            # both rooted
            if self.rootNode != other.rootNode: return False
        elif hasattr(self,"arbitraryNode") and  hasattr(other,"arbitraryNode"):
            # both unrooted
            if self.arbitraryNode != other.arbitraryNode: return False
        else:
            # different types of tree
            return False

        # if we made it here, they are same type of tree, and basis
        # node matches

        # check node connections
        if sorted(self.nodeConnectD.items()) != sorted(other.nodeConnectD.items()):
            return False

        # check branch lengths if they exist
        if self.branchLenD == other.branchLenD:
            return True
        else:
            # don't match
            if self.branchLenD == None or other.branchLen == None:
                # one is None
                return False
            else:
                # both are dicts. Check more carefully
                selfL = []
                for brT,brLen in self.branchLenD.items():
                    selfL.append((brT,round(brLen,4)))
                otherL = []
                for brT,brLen in other.branchLenD.items():
                    otherL.append((brT,round(brLen,4)))

                selfL.sort()
                otherL.sort()
                return selfL == otherL
        
class Rtree(Tree):
    def __init__(self, nodeConnectD=None,rootNode=None):
        '''Initialize a rooted tree object.'''
        super().__init__(nodeConnectD)
        self.rootNode = rootNode

        if self.nodeConnectD != None:
            self.preOrderT = self.__traversePreOrder__(self.rootNode)
        
    def fromNewickFileLoadSpeciesTree(self,treeFN,outGroupTaxaL=None,includeBrLen=False):
        '''Populate attributes based on newick file in treeFN. This method
assumes we are working with a species tree. It must be rooted,
bifurcating and have named internal nodes. An optional argument
outGroupTaxaL can be used to provide a list of outgroups to root the
tree. If outGroupTaxaL is present, we assume we are being asked to
prepare a species tree (e.g. from ASTRAL output) and we do other
preparation such as naming internal nodes).

        '''
        if outGroupTaxaL == None:
            bpTree = Phylo.read(treeFN, 'newick', rooted=True)
        else:
            bpTree = Phylo.read(treeFN, 'newick')
            bpTree = trees.prepareTree(bpTree,outGroupTaxaL)
            
        self.__checkSpeciesTree__(bpTree)

        nodeConnectD,branchLenD = self.__bioPhyloToNodeConnectD__(bpTree)

        if includeBrLen:
            self.populateAttributes(nodeConnectD,bpTree.root.name,branchLenD)
        else:
            self.populateAttributes(nodeConnectD,bpTree.root.name)

    def fromNewickFile(self,treeFN,includeBrLen=False):
        '''Populate attributes based on newick file in treeFN. Assumes a rooted tree.
        '''
        bpTree = Phylo.read(treeFN, 'newick', rooted=True)
        nodeConnectD,branchLenD = self.__bioPhyloToNodeConnectD__(bpTree)
        if includeBrLen:
            self.populateAttributes(nodeConnectD,bpTree.root.name,branchLenD)
        else:
            self.populateAttributes(nodeConnectD,bpTree.root.name)
            
    def toNewickStr(self,includeBrLength=False,nodeLabelD=None):
        '''Output a newick string.'''
        if self.nodeConnectD == None:
            return "()"
        else:
            return self.__traverseForNewickStr__(self.rootNode,ROOT_PARENT_NAME,includeBrLength,nodeLabelD)
        
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
        newConnecT = (ROOT_PARENT_NAME,)+oldConnecT[1:]
        subD[node] = newConnecT

        # make tree
        outTreeO = Rtree()

        # if we have branch lens, include these
        if self.branchLenD == None or len(subD)<2:
            # if self has no branch lengths, or if subtree is a tip, then branchLenD will be None
            outTreeO.populateAttributes(subD,node)
        else:
            subBranchLenD = {}
            for tempNode,connecT in subD.items():
                parentNode = connecT[0]
                if parentNode != ROOT_PARENT_NAME:
                    subBranchLenD[(parentNode,tempNode)] = self.branchLenD[(parentNode,tempNode)]

            outTreeO.populateAttributes(subD,node,subBranchLenD)

        return outTreeO
        
    def createSubtreeD(self):
        '''Get all subtrees and put them in a dict keyed by node name.'''
        subtreeD = {}
        for node in self.preorder():
            subtree = self.subtree(node)
            subtreeD[node] = subtree
        return subtreeD
            
    def children(self,node):
        '''Return tuple of children of node.'''
        connecT = self.nodeConnectD[node]
        return connecT[1:]
        
    def ancestors(self,node):
        '''Return a tuple of nodes ancestral to node.'''

        def traverse(D,node):
            connecT=D[node]
            parent = connecT[0]
            if parent == ROOT_PARENT_NAME: # at root
                return []
            else:
                return [parent]+traverse(D,parent)

        return tuple(traverse(self.nodeConnectD,node))

    def getNearestNeighborL(self,leaf):
        '''Given a leaf, return a list containing the other leaf or
    leaves which are most closely related to leaf.'''
        parent = self.getParent(leaf) # assume leaf is in tree
        subRtreeO = self.subtree(parent)
        leafL = list(subRtreeO.leaves())
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

    def calcDistanceBetweenNodes(self, ancNode, descNode):
        ''' Takes in two nodes in the tree and calculates the distance 
        assuming that the first node is ancestral'''
        if ancNode == descNode:
            return 0
        else:
            currentNode = descNode
            pairsLeadingToAnc = []
            while currentNode != ancNode:
                if currentNode == self.rootNode:
                    return float('inf')
                higherNode = self.getParent(currentNode)
                pairsLeadingToAnc.append((currentNode,higherNode))
                currentNode = higherNode
            distSum = 0
            for pair in pairsLeadingToAnc:
                if pair in self.branchLenD:
                    distSum += self.branchLenD[pair]
                else:
                    distSum += self.branchLenD[(pair[1],pair[0])]
            return distSum

    def prune(self,leafL):
        '''Prune down to a tree with only the leaves in leafL, combining
branches and removing nodes as needed. Returns a new Rtree.
        '''

        # funcs
        def traverseBlD(subTreeO,activeNodeS,blD,pruneParent,previousBrLen,node):
            '''Traverse subtree, keeping only nodes in activeNodeS and adjusting
            branch lengths (add to blD).
            '''
            if subTreeO.isLeaf(node):
                if pruneParent != ROOT_PARENT_NAME:
                    blD[(pruneParent,node)] = previousBrLen + subTreeO.branchLenD[(subRtreeO.getParent(node),node)]
                return blD
            else:
                activeChildrenL = [child for child in subTreeO.children(node) if child in activeNodeS]
                if len(activeChildrenL) > 1:
                    # we don't need to remove this node
                    if pruneParent != ROOT_PARENT_NAME:
                        # here we really want to pass in subRtreeO.getParent(node)
                        # which is the immediate last one, even if its not included in the pruned tree
                        blD[(pruneParent,node)] = previousBrLen + subTreeO.branchLenD[(subRtreeO.getParent(node),node)]
                    # recurse
                    for child in activeChildrenL:
                        blD = traverseBlD(subTreeO,activeNodeS,blD,node,0,child)
                else:
                    # there is only one active child, and we should remove this node
                    child = activeChildrenL[0]
                    brLenToPassAlong = previousBrLen + subTreeO.branchLenD[(subRtreeO.getParent(node),node)]
                    blD = traverseBlD(subTreeO,activeNodeS,blD,pruneParent,brLenToPassAlong,child)
                    
                return blD
            
        def buildNcD(blD,ncD,parent,node):
            '''Fillout ncD (nodeConnectD) based on the branches and lengths in blD.'''
            # get children of node and put in connecT
            childL = []
            for branchPair in blD:
                if branchPair[0] == node:
                    childL.append(branchPair[1])

            connecT = tuple([parent] + childL)
            ncD[node] = connecT

            # recurse. childL will be empty if node a tip.
            for child in childL:
                ncD = buildNcD(blD,ncD,node,child)

            return ncD
                
        # main section
        mrcaNode = self.findMrca(leafL)
        subRtreeO = self.subtree(mrcaNode)

        activeNodeS = set() # leafL + all ancestors
        for leaf in leafL:
            activeNodeS.add(leaf)
            for node in subRtreeO.ancestors(leaf):
                activeNodeS.add(node)

        blD = {}
        blD = traverseBlD(subRtreeO,activeNodeS,blD,subRtreeO.getParent(mrcaNode),0,mrcaNode)
        ncD = {}
        ncD = buildNcD(blD,ncD,ROOT_PARENT_NAME,mrcaNode)
        
        outRtreeO = Rtree()
        outRtreeO.populateAttributes(ncD,mrcaNode,blD)

        return outRtreeO

    def unroot(self):
        '''Return an Utree object corresponding to self. We assume self has at
least two tips.'''
        # funcs
        def traverse(nodeConnectD,branchLenD,parent,node):
            '''Copy entries into nodeConnecteD, starting with node.'''
            if self.isLeaf(node):
                nodeConnectD[node] = (parent,)
            else:
                childT = self.children(node)
                nodeConnectD[node] = (parent,) + childT
                for child in childT:
                    # get branch to child
                    if branchLenD != None:
                        branchLenD[(node,child)] = self.branchLenD[(node,child)]
                    
                    nodeConnectD,branchLenD = traverse(nodeConnectD,branchLenD,node,child)
                
            return nodeConnectD,branchLenD
                
        # unroot main
        nodeConnectD={}
        if self.branchLenD == None:
            branchLenD = None
        else:
            branchLenD = {}

        childT = self.children(self.rootNode)
        if len(childT) < 2:
            raise ValueError("Trying to unroot a single tip tree.")
        elif len(childT) == 2:
            # handle this first branch
            if branchLenD != None:
                branchLenD[(childT[0],childT[1])] = self.branchLenD[(self.rootNode,childT[0])] + self.branchLenD[(self.rootNode,childT[1])]
            # recurse
            nodeConnectD,branchLenD = traverse(nodeConnectD,branchLenD,childT[0],childT[1])
            nodeConnectD,branchLenD = traverse(nodeConnectD,branchLenD,childT[1],childT[0])
            arbitraryNode = childT[0]
        else:
            # Multifurcation at root. Not clear this will ever come up.
            for child in childT:
                nodeConnectD,branchLenD = traverse(nodeConnectD,branchLenD,node,child)
            arbitraryNode = self.rootNode
            
        utreeO = Utree()
        utreeO.populateAttributes(nodeConnectD,arbitraryNode,branchLenD)
        return utreeO
        
    def createDtlorD(self,spTreeB):
        '''Return a dict in the format expected by the dtlor reconciliation
code. Key is a node, value is (parent, node, child1, child2). For now
this only works for bifurcating trees. spTreeB is a boolean that is
True if this is a species tree, and False if it is for a gene
tree.
        '''
        if self.multifurcatingNodes() != []:
            raise ValueError("This tree is not bifurcating.")
        
        dtlorD = OrderedDict()
        for node in self.preorder():
            connecT = self.nodeConnectD[node]
            if len(connecT) == 1:
                # leaf
                parent = connecT[0]
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
    
    def __checkSpeciesTree__(self,bpTree):
        '''Check that a biopython tree is rooted, bifurcating, and has named
    internal nodes. Throw error if not. Returns None.
        '''
        # check its bifurcating. This actually ignores the root, but checks
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

        # Make sure all names unique
        nodeL = [n.name for n in bpTree.get_nonterminals() + bpTree.get_terminals()]
        if len(nodeL) != len(set(nodeL)):
            raise ValueError("All nodes of the input tree (tips and interal) must have unique names.")
        
    def __bioPhyloToNodeConnectD__(self,bpTree):
        '''Convert a biopython tree object to a node connection dict with keys
    that are nodes, and values that are (parent, child1, child2...).'''
        nodeConnectD = {}
        branchLenD = {}
        # parent of root is ROOT_PARENT_NAME
        return self.__bioPhyloCladeToNodeConnectD__(bpTree.root,nodeConnectD,branchLenD,ROOT_PARENT_NAME)

    def __bioPhyloCladeToNodeConnectD__(self,clade,nodeConnectD,branchLenD,parent):
        '''Recursive helper to convert a biopython clade object.
        '''
        # get branch, but only if its not branch leading to root
        if parent != ROOT_PARENT_NAME: branchLenD[(parent,clade.name)] = clade.branch_length

        # fill nodeConnectD
        if clade.is_terminal():
            nodeConnectD[clade.name] = (parent,)
            return nodeConnectD,branchLenD
        else:
            nodeConnectD[clade.name] = (parent,clade[0].name,clade[1].name)
            for iterClade in clade:
                nodeConnectD,branchLenD = self.__bioPhyloCladeToNodeConnectD__(iterClade,nodeConnectD,branchLenD,clade.name)
            return nodeConnectD,branchLenD

    def __repr__(self):
        return "Rtree: "+self.toNewickStr()


class Utree(Tree):
    def __init__(self, nodeConnectD=None):
        '''Initialize an unrooted tree object.'''
        super().__init__(nodeConnectD)
        
        if self.nodeConnectD != None:
            # arbitraryNode provides a starting point for some methods
            self.arbitraryNode = self.internals()[0] if len(self.internals())>0 else self.leaves()[0]
            self.preOrderT = self.__traversePreOrder__(self.arbitraryNode)
            self.branchPairT = self.__createBranchPairT__()
            
    def fromNewickFile(self,treeFN):
        '''Populate attributes based on newick file in treeFN. This method
assumes we are working with an unrooted gene tree that does not have
named interal nodes (we will create those).

        '''
        bpTree = Phylo.read(treeFN, 'newick', rooted=False)
        # handle special case of one tip tree
        if bpTree.count_terminals()==1:
            # special case, one tip tree
            nodeConnectD={bpTree.get_terminals()[0].name:()}
            branchLenL=[('','',None)]
        else:
            iNodeNum,nodeConnectD,branchLenL = self.__bioPhyloToNodeConnectD__(bpTree.clade,ROOT_PARENT_NAME,{},0)

        # get arbitrary node
        tempLeavesT,tempInternalsT = self.__updateSecondaryAttributesHelper__(nodeConnectD)
        arbitraryNode = tempInternalsT[0] if len(tempInternalsT)>0 else tempLeavesT[0]

        # branchLenD
        branchLenD = {}
        for parent,child,brLen in branchLenL:
            if parent != ROOT_PARENT_NAME:
                branchLenD[(parent,child)] = brLen
        if all(brLen==None for brLen in branchLenD.values()):
            # if all branches None, don't bother with dict
            branchLenD = None
        
        # put in these attributes
        self.populateAttributes(nodeConnectD,arbitraryNode,branchLenD)

    def toNewickStr(self,includeBrLength=False,nodeLabelD=None):
        '''Output a newick string.'''
        if self.nodeConnectD == None:
            # empty
            return "()"
        elif self.nodeCount() == 1:
            return "("+self.leafNodeT[0]+")"
        elif self.nodeCount() == 2:
            # two tip tree is a special case            
            return "("+self.leafNodeT[0]+","+self.leafNodeT[1]+")"
        else:
            return self.__traverseForNewickStr__(self.arbitraryNode,ROOT_PARENT_NAME,includeBrLength,nodeLabelD)

    def root(self,branchToRootPair):
        '''Root using the tuple branchToRootPair and return an Rtree object.'''

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
        self.__splitNodeConnectD__(self.nodeConnectD,newD,branchToRootPair[0],branchToRootPair[1])
        self.__splitNodeConnectD__(self.nodeConnectD,newD,branchToRootPair[1],branchToRootPair[0])

        # must add the root
        rootNode = "root"
        newD[rootNode] = (ROOT_PARENT_NAME,branchToRootPair[0],branchToRootPair[1])

        # adjust children of root to say root is parent
        newD = updateChildOfRoot(newD,branchToRootPair[0],branchToRootPair[1])
        newD = updateChildOfRoot(newD,branchToRootPair[1],branchToRootPair[0])
        
        return Rtree(newD,rootNode)

    def rootIncludeBranchLen(self,branchToRootPair):
        '''Root using the tuple branchToRootPair and return an Rtree object that
includes branch lengths.'''

        # output tree
        rtreeO = self.root(branchToRootPair)

        # get all branch lengths besides those involving the root
        branchLenD = {}
        for utreeBranchPair in self.branchPairT:

            if utreeBranchPair == branchToRootPair or utreeBranchPair == (branchToRootPair[1],branchToRootPair[0]):
                # this is the branch to be rooted, skip
                continue
            else:
                # make sure its in preorder as defined for rtreeO
                rtreeBranchPair = tuple((nd for nd in rtreeO.preorder() if nd in utreeBranchPair))

                branchLenD[rtreeBranchPair] = self.branchLenD[utreeBranchPair]

        # now get branch lens involving root. Arbitrarily split 50:50.
        utreeBranchToRootLen = self.branchLenD[branchToRootPair]
        branchLenD[("root",branchToRootPair[0])] = 0.5 * utreeBranchToRootLen
        branchLenD[("root",branchToRootPair[1])] = 0.5 * utreeBranchToRootLen
        
        # add to tree
        rtreeO.branchLenD = branchLenD

        return rtreeO
        
    def iterAllRootedTrees(self):
        '''Iterator yielding all possible rooted trees from this unrooted tree.'''
        for branchPair in self.branchPairT:
            yield self.root(branchPair)

    def iterAllRootedTreesIncludeBranchLen(self):
        '''Iterator yielding all possible rooted trees from this unrooted tree.'''
        for branchPair in self.branchPairT:
            yield self.rootIncludeBranchLen(branchPair)

    def split(self,branchPair):
        '''Split on the branch specified by branchPair into two new Utree
objects.

        '''
        
        def subUtree(oldNodeConnectD,oldBranchPairT,oldBranchLenD,node,parentNode):
            '''Given an oldNodeConnectD from a Utree object, and a node and it's
parent, create a new Utree object which is a sub tree. Assumes node is
not a tip.

            '''
            # get subset of oldNodeConnectD, put in D
            D = {}
            self.__splitNodeConnectD__(oldNodeConnectD,D,node,parentNode)

            newBranchLenD = {}
            if len(D[node])>3:
                # it's a multifurcating node, node should remain in D,
                # but we must remove the connection to parentNode
                connecT = D[node]
                newConnecL = []
                for tempNode in connecT:
                    if tempNode != parentNode:
                        newConnecL.append(tempNode)
                D[node] = tuple(newConnecL)

                # create output tree
                newUtreeO = Utree(D)
                # get branch lengths
                for nbp in newUtreeO.branchPairT:
                    # might need to reverse order of nbp to find in oldBranchLenD
                    if nbp in oldBranchLenD:
                        newBranchLenD[nbp] = oldBranchLenD[nbp]
                    else:
                        newBranchLenD[nbp] = oldBranchLenD[(nbp[1],nbp[0])]
                   
            else:
                # bifurcating node
                # get the two branches that connect to node so we can
                # remove node itself.
                _,child1,child2 = D[node]

                assert(_==parentNode) # temp

                # adjust entries for children
                updateChild(D,node,child1,child2)
                updateChild(D,node,child2,child1)
                del D[node]
            
                # create output tree
                newUtreeO = Utree(D)

                # get branch lengths
                for nbp in newUtreeO.branchPairT:
                    if child1 in nbp and child2 in nbp:
                        # this branch len must be made by adding two old ones
                        oldBP1 = [oldBP for oldBP in oldBranchPairT if (node in oldBP and child1 in oldBP)][0]
                        oldBP2 = [oldBP for oldBP in oldBranchPairT if (node in oldBP and child2 in oldBP)][0]
                        brLen = oldBranchLenD[oldBP1] + oldBranchLenD[oldBP2]
                        newBranchLenD[nbp] = brLen
                    else:
                        # might need to reverse order of nbp to find in oldBranchLenD
                        if nbp in oldBranchLenD:
                            newBranchLenD[nbp] = oldBranchLenD[nbp]
                        else:
                            newBranchLenD[nbp] = oldBranchLenD[(nbp[1],nbp[0])]
                        
            newUtreeO.branchLenD = newBranchLenD
            return newUtreeO

        def updateChild(D,node,child1,child2):
            '''Update entries in D for child1 to connect to child2 rather than
node.'''
            connecL = list(D[child1])
            connecL[connecL.index(node)] = child2
            D[child1] = tuple(connecL)
            return

        # main part of split
        if self.isLeaf(branchPair[0]) and self.isLeaf(branchPair[1]):
            return Utree({branchPair[0]:()}),Utree({branchPair[1]:()}) # branchLenD None by default
        elif self.isLeaf(branchPair[0]):
            # one side is leaf
            restUtreeO = subUtree(self.nodeConnectD,self.branchPairT,self.branchLenD,branchPair[1],branchPair[0])
            return Utree({branchPair[0]:()}),restUtreeO
        elif self.isLeaf(branchPair[1]):
            restUtreeO = subUtree(self.nodeConnectD,self.branchPairT,self.branchLenD,branchPair[0],branchPair[1])
            return restUtreeO,Utree({branchPair[1]:()})
        else:
            # internal
            aUtreeO = subUtree(self.nodeConnectD,self.branchPairT,self.branchLenD,branchPair[1],branchPair[0])
            bUtreeO = subUtree(self.nodeConnectD,self.branchPairT,self.branchLenD,branchPair[0],branchPair[1])
            return aUtreeO,bUtreeO

    def makeDistanceMatrix(self) -> dict:
        """ makeDistanceMatrix:
                Accepts a Utree object as input. Constructs and returns a matrix 
                of the cophenetic distances for all of the pairwise combinations of
                the tips. Returns the matrix as a dictionary keyed by tuples conta-
                ining each pairwise combination.
        """
        # get all tips as a list
        tipsL = list(self.leafNodeT)
        # build a dictionary of all pairwise combinations
        tableD = dict()
        for tip1 in tipsL:
            for tip2 in tipsL:
                # self-vs-self comparisons are on the diagonal; they should be zero
                if tip1 == tip2:
                    distance = 0.0
                # calculate the cophenetic distances for each comparison
                else:
                    # calculate the distance between the two tips
                    path = self.__treePathfinder__(tip1, tip2)
                    distance = self.__sumBranchLengths__(path)
                # save the calculated distance in the tableD
                tableD[(tip1, tip2)] = distance
        return tableD
    
    def __treePathfinder__(self, origin:str, end:str, previous:str=None) -> list:
        """ __treePathfinder__:
                Accepts a Utree object, a string indicating the starting node/tip,
                a string indicating the ending node/tip, and a string indicating 
                the previous origin (used for recursion) as inputs. Recurses throu-
                gh the the tree to find the branches (indicated by pairs of nodes)
                that will connect the origin to the end within the tree. Returns a 
                list of tuples where each tuple contains a pair of nodes.
        """
        # get all the nodes/tips that are connected to the origin
        allConnxns = self.nodeConnectD[origin]
         # base case: the end is connected to the origin
        if end in allConnxns:
            return [(origin, end)]
        # iterate through the origin's connections
        outL = list()
        for connxn in allConnxns:
            # only process connections that aren't leaves or the previous node
            if not self.isLeaf(connxn) and connxn != previous:
                # recurse on connection towards the specified 'end'
                path = self.__treePathfinder__(connxn, end, previous=origin)
                # if a path was found ...
                if path != []:
                    # ... then append the origin and its connection to the list ...
                    outL.append((origin, connxn))
                    # ... and add the path that was found to the list
                    outL += path
        # returns an empty list if nothing was found or if nothing to do
        return outL

    def __sumBranchLengths__(self, pathL:list) -> float:
        """ __sumBranchLengths__:
                Accepts a Utree object and a list of node-pairs (tuples) as inputs.
                Calculates and returns the total distance represented by the list
                of node-pairs.
        """
        # initialize the sum
        totalLen = 0.0
        # for each node-pair in the list
        for nodePair in pathL:
            # flip the order if the pair is not a key in branchLenD
            if nodePair not in self.branchLenD.keys():
                nodePair = (nodePair[1], nodePair[0])
            # add the length to the running total
            totalLen += self.branchLenD[nodePair]
        return totalLen

    def __bioPhyloToNodeConnectD__(self,bpClade,parentNodeStr,nodeConnectD,iNodeNum):
        '''Convert a biopython clade object to a node connection dict with
keys that are nodes, and values that are the nodes connected to. For
internals (node1, node2, node3...). We assume bpClade does not have
named internal nodes, and we name them here: g0, g1 etc.

        '''

        if bpClade.is_terminal():
            # parentNodeStr will never be ROOT_PARENT_NAME in this case
            nodeConnectD[bpClade.name] = (parentNodeStr,)
            branchLenL = [(parentNodeStr,bpClade.name,bpClade.branch_length)]

        elif parentNodeStr == ROOT_PARENT_NAME and len(bpClade) == 2:
            # biopython has put 2 branches only at the base. We will
            # operate over this branch, joining the two nodes with no
            # node in between.

            # get branch length for new branch with root removed
            if bpClade[0].branch_length == None or bpClade[1].branch_length == None:
                brLen = None
            else:
                brLen = bpClade[0].branch_length + bpClade[1].branch_length
            
            # handle cases where one or two are terminal
            if bpClade[0].is_terminal() and bpClade[1].is_terminal():
                # both terminal
                iNodeNum,nodeConnectD,_ = self.__bioPhyloToNodeConnectD__(bpClade[0],bpClade[1].name,nodeConnectD,iNodeNum)
                iNodeNum,nodeConnectD,_ = self.__bioPhyloToNodeConnectD__(bpClade[1],bpClade[0].name,nodeConnectD,iNodeNum)
                
                branchLenL = [(bpClade[0].name,bpClade[1].name,brLen)]
                
            elif bpClade[0].is_terminal():
                # 0 is terminal, 1 not
                # parent node for 0 will be based on current iNodeNum
                iNodeNum,nodeConnectD,_ = self.__bioPhyloToNodeConnectD__(bpClade[0],"g"+str(iNodeNum),nodeConnectD,iNodeNum)
                iNodeNum,nodeConnectD,branchLenL = self.__bioPhyloToNodeConnectD__(bpClade[1],bpClade[0].name,nodeConnectD,iNodeNum)

                # fix len of br made when root removed
                for i in range(len(branchLenL)):
                    if bpClade[0].name in branchLenL[i]:
                        # only be 1 branch w/this in it. Update the
                        # len to proper one calculated above
                        branchLenL[i] = (branchLenL[i][0],branchLenL[i][1],brLen)
                
            elif bpClade[1].is_terminal():
                # 1 is terminal, 0 not
                # parent node for 1 will be based on current iNodeNum
                iNodeNum,nodeConnectD,_ = self.__bioPhyloToNodeConnectD__(bpClade[1],"g"+str(iNodeNum),nodeConnectD,iNodeNum)
                iNodeNum,nodeConnectD,branchLenL = self.__bioPhyloToNodeConnectD__(bpClade[0],bpClade[1].name,nodeConnectD,iNodeNum)
                # fix len of br made when root removed
                for i in range(len(branchLenL)):
                    if bpClade[1].name in branchLenL[i]:
                        branchLenL[i] = (branchLenL[i][0],branchLenL[i][1],brLen)

            else:
                # neither is terminal

                # get the iNodeNum for the first internal node in 0
                iNodeNumAtBaseOf0 = iNodeNum
                
                # figure out what the parent for 0 will be by making a
                # dummy call. iNodeNum afterward will be the first
                # internal node in 1.
                iNodeNumAtBaseOf1,_,_ = self.__bioPhyloToNodeConnectD__(bpClade[0],"dummy",{},iNodeNum)

                # now make real calls
                branchLenL = []
                iNodeNum,nodeConnectD,tempBranchLenL = self.__bioPhyloToNodeConnectD__(bpClade[0],"g"+str(iNodeNumAtBaseOf1),nodeConnectD,iNodeNum)
                branchLenL.extend(tempBranchLenL)

                iNodeNum,nodeConnectD,tempBranchLenL = self.__bioPhyloToNodeConnectD__(bpClade[1],"g"+str(iNodeNumAtBaseOf0),nodeConnectD,iNodeNum)

                # there is a redundant entry in tempBranchLenL. Don't
                # include in branchLenL. The entry that has
                # ("g"+str(iNodeNumAtBaseOf0),"g"+str(iNodeNumAtBaseOf1),brLen)
                for brLenT in tempBranchLenL:
                    if brLenT[0] == "g"+str(iNodeNumAtBaseOf0) and brLenT[1] == "g"+str(iNodeNumAtBaseOf1):
                        pass
                    elif brLenT[1] == "g"+str(iNodeNumAtBaseOf0) and brLenT[0] == "g"+str(iNodeNumAtBaseOf1):
                        # keep this, updating len of br
                        newBrLenT = (brLenT[1],brLenT[0],brLen) # brLen is correct len of br where root removed
                        branchLenL.append(newBrLenT)
                    else:
                        branchLenL.append(brLenT)
        else:
            # general case
            connecL = []
            branchLenL = []
            if parentNodeStr != ROOT_PARENT_NAME:
                connecL.append(parentNodeStr)
            # if we're just starting, it's ROOT_PARENT_NAME, and we
            # don't want to put that in the list of connections
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
                iNodeNum,nodeConnectD,tempBranchLenL = self.__bioPhyloToNodeConnectD__(childClade,thisNodeStr,nodeConnectD,iNodeNum)
                branchLenL.extend(tempBranchLenL)
                
            # add in this connection for this node
            nodeConnectD[thisNodeStr] = tuple(connecL)

            # get branch len
            if hasattr(bpClade,'branch_length'):
                brLen = bpClade.branch_length
            else:
                brLen = None
            branchLenL.append((parentNodeStr,thisNodeStr,brLen))

            
        return iNodeNum,nodeConnectD,branchLenL
    
    def __splitNodeConnectD__(self,D,newD,node,parentNode):
        '''Get nodeConnectD for the part of the tree defined by node, and in
the oposite direction from parentNode. Put in newD. Returns None.

        '''
        oldConnecT=D[node]

        # make sure parent node is first
        assert(parentNode in oldConnecT)
        connecL = [parentNode]
        for tempNode in oldConnecT:
            if tempNode != parentNode:
                connecL.append(tempNode)
        connecT = tuple(connecL)

        newD[node] = connecT # store

        if len(connecT)==1:
            return
        else:
            for child in connecT:
                if child != parentNode:
                    self.__splitNodeConnectD__(D,newD,child,node)
            return
        
    def __repr__(self):
        return "Utree: "+self.toNewickStr()
