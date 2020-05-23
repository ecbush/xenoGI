# functions for loading, manipulating and creating phylogenetic trees
from Bio import Phylo
import sys, os, glob, subprocess, shutil
from multiprocessing import Pool
from xenoGI import xenoGI, blast, scores, parameters, fasta, genomes

#### General tree related functions

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
    
def prune(tree,strainsToKeep):
    '''Given a tree and a group of strains strainsToKeep (which we assume
are a subset of the strains in tree), prune to yield the tree with
only those strains. We will assume tree has no branch lengths, and
will thus always put None in the 4 position of the pruned tree.
    '''
    if tree[1]==():
        if tree[0] in strainsToKeep:
            return (tree[0],(),(),None)
        else:
            return None
    else:
        l = prune(tree[1],strainsToKeep)
        r = prune(tree[2],strainsToKeep)

        if l != None and r != None:
            return (tree[0],l,r,None)
        elif l == None:
            return r
        else:
            return l
        
def isRootNode(tree,mrcaNum):
    '''Is mrcaNum the root node?'''
    # root node is tree[0] in our tuple trees
    return mrcaNum == tree[0]
    
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
    tupleTree=bioPhyloToTupleTree(bpTree)
    return tupleTree

def createSubtreeL(tree):
    '''Return a list containing all subtrees.'''
    if tree[1]==():
        return [tree]
    else:
        return [tree]+createSubtreeL(tree[1]) + createSubtreeL(tree[2])

def createSubtreeD(tree):
    '''Get all subtrees and put them in a dict keyed by root node name.'''
    D = {}
    for subtree in createSubtreeL(tree):
        D[subtree[0]] = subtree
    return D
    
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

def getRootFocalCladeFromOutgroup(tree,outGroup):
    '''Given a rooted tree and an outgroup, find the
    rootFocalClade. Assumes the outgroup is present in tree.
    '''
    left = tree[1]
    right = tree[2]
    if isSpeciesPresent(right,outGroup):
        return left[0]    
    elif isSpeciesPresent(left,outGroup):
        return right[0]
    else:
        raise ValueError("Error trying to find rootFocalClade, outGroup not present in provided tree.")
    
def isSpeciesPresent(tree,species):
    '''Return True is species in tree, False otherwise.'''
    if tree[1]==():
        if tree[0]==species: return True
        else: return False
    else:
        l=isSpeciesPresent(tree[1],species)
        r=isSpeciesPresent(tree[2],species)
        return l or r
        
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

#### Functions for creating gene and species trees

def makeSpeciesTree(paramD,aabrhHardCoreL,genesO):
    '''Create a species tree based on the orhtolog sets in aabrhHardCoreL.'''

    ## set up
    # create work dir if it doesn't already exist
    workDir = paramD['makeSpeciesTreeWorkingDir']
    if glob.glob(workDir)==[]:
        os.mkdir(workDir)

    numProcesses = paramD['numProcesses']
    javaPath = paramD['javaPath']
    astralPath = paramD['astralPath']
    astralTreeFN = paramD['astralTreeFN']
    gtFileStem = 'aabrhHardCoreFam'
    allGtFilePath = os.path.join(workDir,gtFileStem+'*.tre')
    aabrhHardCoreGeneTreesFN = paramD['aabrhHardCoreGeneTreesFN']
    outSpeciesTreeFN = paramD['treeFN'] # for main output
    deleteSpeciesTreeWorkingDir = paramD['deleteSpeciesTreeWorkingDir']
    outGroupTaxaL = [paramD['outGroup']]

    # if tree file already exists, throw error
    if os.path.isfile(paramD['treeFN']):
        raise IOError("The tree file " + paramD['treeFN'] + " already exists.")

    ## make gene tree for each aabrh hard Core set
    # add numbering to list
    newAabrhHardCoreL = []
    for orthoNum,orthoT in enumerate(aabrhHardCoreL):
        newAabrhHardCoreL.append((orthoNum,orthoT))
    
    makeGeneTrees(paramD,True,genesO,workDir,gtFileStem,newAabrhHardCoreL)

    ## run Astral on gene trees

    # concatenate all gene tree files
    with open(aabrhHardCoreGeneTreesFN, 'w') as aabrhHardCoreGeneTreesF:
        for fn in sorted(glob.glob(allGtFilePath)):
            with open(fn) as infile:
                aabrhHardCoreGeneTreesF.write(infile.read())
    
    subprocess.call([javaPath, '-jar',astralPath, '-i', aabrhHardCoreGeneTreesFN, '-o',astralTreeFN],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    ## load and root the tree.
    bpTree = Phylo.read(astralTreeFN, 'newick')
    tupleTree=prepareTree(bpTree,outGroupTaxaL)
    writeTreeNoBrLen(tupleTree,outSpeciesTreeFN) # final output to 
    
    ## delete the working directory
    if deleteSpeciesTreeWorkingDir:
        shutil.rmtree(workDir)

def makeGeneTrees(paramD,strainHeader,genesO,workDir,gtFileStem,orthoTL):
    '''Given a list of ortho groups, orthoTL, make a gene tree for each
and put in workDir. strainHeader is a boolean that if True means we
put the strain in the header for alignments (so we'll end up with that
for the names of the tips of the tree). gtFileStem gives the stem of
the name for the gene tree files. orthoTL is a list of
(orthoGroupNum,orthoT) where orthoT has the genes in an ortholog
group.

    '''
    numProcesses = paramD['numProcesses']    
    musclePath = paramD['musclePath']
    fastTreePath = paramD['fastTreePath']
    
    ## make gene trees
    argumentL = [([],[],paramD,strainHeader,genesO,workDir,gtFileStem,musclePath,fastTreePath) for i in range(numProcesses)]

    # set up argumentL
    counter = 0
    for orthoGroupNum,orthoT in orthoTL:
        # we use counter for placememnt in argumentL since at least in
        # the case of single gene families, some will be missing
        orthoGroupNumStr = str(orthoGroupNum).zfill(6) # pad w/ 0's so ls will display in right order
        argumentL[counter%numProcesses][0].append(orthoGroupNumStr)
        argumentL[counter%numProcesses][1].append(orthoT)
        counter += 1
        
    # run
    with Pool(processes=numProcesses) as p:
        for _ in p.imap_unordered(makeOneGeneTreeGroup, argumentL):
            pass

    return
        
def makeOneGeneTreeGroup(argT):
    '''Wrapper for multiprocessing. Given a group of orthoT's to work on,
calls makeOneGeneTree on each.'''
    
    orthoGroupNumStrL,orthoTL,paramD,strainHeader,genesO,workDir,gtFileStem,musclePath,fastTreePath = argT

    # get set of all genes in orthoTL to restrict size of seqD's
    orthoGenesS=set()
    for orthoT in orthoTL:
        orthoGenesS.update(orthoT)
    
    # load protein and dna sequences
    protSeqD=genomes.loadSeq(paramD, '_prot.fa',genesS=orthoGenesS)

    if paramD['dnaBasedGeneTrees'] == True:
        dnaSeqD = genomes.loadSeq(paramD, '_dna.fa',genesS=orthoGenesS)
    else:
        dnaSeqD = {}

    # make trees
    for i in range(len(orthoTL)):
        orthoGroupNumStr = orthoGroupNumStrL[i]
        orthoT = orthoTL[i]
        makeOneGeneTree(orthoGroupNumStr,orthoT,strainHeader,genesO,protSeqD,dnaSeqD,workDir,gtFileStem,musclePath,fastTreePath)
    
def makeOneGeneTree(orthoGroupNumStr,orthoT,strainHeader,genesO,protSeqD,dnaSeqD,workDir,gtFileStem,musclePath,fastTreePath):
    ''' Makes one gene tree from an ortho list. If dnaSeqD is empty, uses protein only.'''
 
    # get temp align file names
    inTempProtFN=os.path.join(workDir,"tempProt"+orthoGroupNumStr+".fa")
    outAlignFN=os.path.join(workDir,"align"+orthoGroupNumStr+".afa")

    # align
    alignOneOrthoT(orthoT,strainHeader,musclePath,inTempProtFN,outAlignFN,protSeqD,dnaSeqD,genesO)
    os.remove(inTempProtFN) # remove unaligned prot file
    
    # make gene tree
    geneTreeFN = os.path.join(workDir,gtFileStem+orthoGroupNumStr+".tre")
    if dnaSeqD == {}:
        # using protein
        subprocess.call([fastTreePath, '-out',geneTreeFN,outAlignFN],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    else:
        # using dna
        subprocess.call([fastTreePath,'-gtr','-nt','-out',geneTreeFN,outAlignFN],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

def alignOneOrthoT(orthoT,strainHeader,musclePath,inProtFN,outAlignFN,protSeqD,dnaSeqD,genesO):
    '''Given genes in a single ortholog set, align them with muscle. If
dnaSeqD is empty, uses protein only.'''

    # write prots we want to align to temp file
    writeFasta(inProtFN,orthoT,strainHeader,genesO,protSeqD)

    # align proteins
    retCode = subprocess.call([musclePath, '-in' ,inProtFN, '-out', outAlignFN],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    if retCode != 0:
        raise OSError("Alignment failed for "+inProtFN)
    
    if dnaSeqD != {}:
        # back align to get dna alignment, overwriting protein alignment.
        protAlignL = []
        for header,alignedProtSeq in fasta.load(outAlignFN):
            if strainHeader:
                protAlignL.append((int(header.rstrip().split()[1]),header,alignedProtSeq))
            else:
                protAlignL.append((int(header.rstrip()[1:]),header,alignedProtSeq))
        backAlign(outAlignFN,protAlignL,dnaSeqD,genesO)

    return
    
def writeFasta(inProtFN,orthoT,strainHeader,genesO,seqD):
    '''Writes a fasta block. seqs are specified in orthoT, and obtained
from seqD. If strainHeader is True, then we put strain whitespace gene
number in header. (In this case, the tree will end up with strain name
on its tips). Otherwise, only gene number.

    '''
    with open(inProtFN,"w") as f:
        for geneNum in orthoT:
            if strainHeader:
                strainName = genesO.numToStrainName(geneNum)
                f.write(">" + strainName + ' '  + str(geneNum) + "\n")
            else:
                f.write(">" + str(geneNum) + "\n")

            f.write(seqD[geneNum] + "\n")
            f.write("\n")

def backAlign(outAlignFN,protAlignL,dnaSeqD,genesO):
    '''Prints the nucleotide alignments given by one block of protein
    alignments.'''
    with open(outAlignFN,'w') as outf:
        printBlock=''
        for geneNum,header,protSeq in protAlignL:
            dnaSeq = dnaSeqD[geneNum]
            lenProtein = len(protSeq) - protSeq.count('-')
            if (lenProtein + 1) * 3 != len(dnaSeq):
                # dna has stop codon                          
                raise IndexError("Lengths of dna and protein do not correspond.")

            strainName = genesO.numToStrainName(geneNum)
            # add header
            printBlock += header + '\n'
            # adding gaps to the nucleotide sequence corresponding to gaps in 
            # the protein sequence and printing the nucleotide sequence
            printBlock += fixSeq(dnaSeq, protSeq) + "\n"
        print(printBlock,file=outf)

def fixSeq(sequence, protAlignment):
    """takes a nucleotide sequence and adds in gaps based on the corresponding
    protein alignment."""
    result = ""
    x = 0
    y = 0

    # while there are still nucleotides in the sequence
    while y != len(protAlignment):

        # if there is a gap in the protein alignment, add gaps to the nucleotide
        # sequence
        if protAlignment[y] == '-':
            result += "---"
            y += 1

        # if there is an amino acid, use the next codon
        else:
            result += sequence[3*x:3*(x+1)].upper()
            x += 1
            y += 1
    return result

## Functions to root, strip branch lengths and name internal nodes

def prepareTree(bpTree,outGroupTaxaL):
    '''Given a biopython unrooted tree object (e.g. ASTRAL output) root,
name the internal nodes, and strip branch lengths (if
any). outGroupTaxaL provides the outgroups for rooting. Return
resulting tree in four tuple form.

    '''
    if not bpTree.is_bifurcating():
        raise ValueError("This tree is not bifurcating. A trifurcation at the root is permitted (to signify it's unrooted) and is not the problem here. There cannot be other trifurcations though, and there are here.")

    # strip confidence values in internal nodes (which mess up our
    # names)
    for n in bpTree.get_nonterminals():
        n.confidence = None
    
    bpTree = rootTree(bpTree,outGroupTaxaL)
    bpTree = nameInternalNodes(bpTree)

    checkTree(bpTree)

    tupleTree=bioPhyloToTupleTree(bpTree)
     

    return tupleTree
    
def rootTree(tree,outGroupTaxaL):
    '''Given a bio python tree, root it. outGroupTaxaL should consist of
all the species in the desired outgroup clade.'''

    tree.root_with_outgroup(*outGroupTaxaL)
    # rooting with outgroups with multiple species seems to leave a
    # trifurcating root.
    if len(tree.root)==2:
        return tree
    else:
        return makeRootBifurcation(tree,outGroupTaxaL)        

def makeRootBifurcation(tree,outGroupTaxaL):
    '''Given a tree that is not bifurcating at root, if it is trifurcating convert to bifurcating. Otherwise thow error.'''

    if len(tree.root)!=3:
        raise ValueError("This tree has a strange number of branches at the root (not 2 or 3).")

    # find the two in outGroupTaxaL, and the one not in it
    inL=[]
    notIn=None
    for cld in tree.root:
        if tipsInOutgroupL(cld,outGroupTaxaL):
            inL.append(cld)
        else:
            notIn = cld

    if len(inL)!=2:
        raise ValueError("There should be two clades of this trifurcating tree which are in outGroupTaxaL, but there are not.")

    # We should take the branch length of notIn and split it up, one
    # half going to notIn, the other half to the new inL clade

    brToSplit=notIn.branch_length
    notIn.branch_length = brToSplit/2
    
    c1=Phylo.Newick.Clade(clades=inL,branch_length=brToSplit/2)
    c2=Phylo.Newick.Clade(clades=[c1,notIn],branch_length=0)
    outTree = Phylo.Newick.Tree(root=c2,rooted=True)

    return outTree
    
def tipsInOutgroupL(clade,outGroupTaxaL):
    '''Return True if all of the tips in clade are in outGroupTaxaL False
otherwise..'''
    for tipName in cladeTipNames(clade):
        if tipName not in outGroupTaxaL:
            return False
    return True
    
def cladeTipNames(clade):
    '''Get the names of all the tips in this clade.'''
    L=[]
    for tipC in clade.get_terminals():
        L.append(tipC.name)
    return L

def nameInternalNodes(tree):
    '''Given a bio python tree without named internal nodes, name
them. This obliterates anything names that might have been inside
already (e.g. if confidence values have been put in the name field).'''
    for j,n in enumerate(tree.get_nonterminals()):
        n.name="i"+str(j)
    return tree 

def tupleTree2NoBrLenNewick(tree):
    '''Convert a four tuple based tree (root,left,right,branchLen) into a
newick formated string, without branch lengths.'''
    if tree[1]==():
        return str(tree[0])
    else:
        leftStr=tupleTree2NoBrLenNewick(tree[1])
        rightStr=tupleTree2NoBrLenNewick(tree[2])
        return "("+leftStr+","+rightStr+")"+str(tree[0])

def writeTreeNoBrLen(tree,fileName):
    '''Write tree to fileName (in newick format).'''
    f=open(fileName,"w")
    f.write(tupleTree2NoBrLenNewick(tree)+";")
    f.close()
    
def stripBranchLen(tree):
    '''Remove branch lengths from tree. Modifies in place.'''
    for n in tree.get_nonterminals() + tree.get_terminals():
        del n.branch_length
    return tree

def makeGeneFamilyTrees(paramD,genesO,familiesO,gtFileStem = 'fam'):
    '''Given a families object, create a gene tree for each family.'''

    # create work dir if it doesn't already exist
    workDir = paramD['geneFamilyTreesDir']
    if glob.glob(workDir)==[]:
        os.mkdir(workDir)

    # create gene tuple for each family
    orthoTL = []
    for familyO in familiesO.iterFamilies():
        familyT = tuple(familyO.iterGenes())
        if len(familyT) > 1:
            # don't bother with single gene fams
            orthoTL.append((familyO.famNum,familyT))
            
    # make gene trees
    makeGeneTrees(paramD,False,genesO,workDir,gtFileStem,orthoTL)

    # remove alignments
    for fn in glob.glob(os.path.join(workDir,"align*.afa")):
        os.remove(fn)
    
    return

#### Gene tree related functions

def loadOneGeneTree(treeFN):
    '''Load an unrooded gene tree. If it is bifurcating arbitrarily root,
convert to tuple tree format and return. If it is multifurcating
return None.
    '''

    bpTree = Phylo.read(treeFN, 'newick', rooted=False)
    bpTree = nameInternalNodes(bpTree)

    # root arbitrarily
    try:
        bpTree.root_at_midpoint()
    except:
        bpTree.root_with_outgroup(bpTree.get_terminals()[0])
        
    # name the root
    bpTree.root.name = "root"

    # if tips have numbers for names, convert to int
    for clade in bpTree.get_terminals():
        if clade.name.isdigit():
            clade.name = int(clade.name)
    
    # if multifurcating return None
    if not bioPhyloIsBinary(bpTree):
        return None
    else:
        # otherwise convert to tupleTree
        tree = bioPhyloToTupleTree(bpTree)
        return tree

def bioPhyloIsBinary(bpTree):
    """
    Checks if a rooted biopython tree is binary.
    """
    # Had some problems with bpTree.is_bifurcating, so wrote this
    preorder=list(bpTree.find_clades(order='preorder'))
    for node in preorder:
        children=[child.name for child in node]
        if len(children)==0 or len(children)==2:
            continue
        else:
            return False
    return True

#### Getting all rootings

def createLocusMapForRootingD(geneTree,locusMapForRootingD):
    """
    Try assigning syntenic locations to internal nodes 
    in the input (arbitrarily rooted) tree. If a loc can be assigned, we are not
    going to root on branches where both terminals are assigned loc
    """
    if geneTree[0] in locusMapForRootingD:
        return locusMapForRootingD
    else:
        # always an internal node since locusMapForRootingD starts w/tips
        # recurse on left and right
        locusMapForRootingD = createLocusMapForRootingD(geneTree[1],locusMapForRootingD)
        locusMapForRootingD = createLocusMapForRootingD(geneTree[2],locusMapForRootingD)
        # handle current node
        if locusMapForRootingD[geneTree[1][0]] == locusMapForRootingD[geneTree[2][0]]:
            locusMapForRootingD[geneTree[0]] = locusMapForRootingD[geneTree[1][0]]
        else:
            locusMapForRootingD[geneTree[0]] = "*"

        return locusMapForRootingD

def left_subtree(tree):
    return tree[1]
def right_subtree(tree):
    return tree[2]
def is_leaf(tree):
    return tree[1]==()
def get_name(tree):
    return tree[0]

def validRooting(tree, locusMapForRootingD):
    if locusMapForRootingD[get_name(tree[1])]==locusMapForRootingD[get_name(tree[2])] and locusMapForRootingD[get_name(tree[1])]!="*":  
        #don't root at a branch when both terminals are imputed the same loc
       return False
    else: return True

def get_all_rerootings(tree, locusMapForRootingD):
    """
    Compute all rerootings of a tree

    Parameters
    --------------------
        tree                   -- a four tuple tree
        locusMapForRootingD   -- a dictionary mapping the gene node name
                                 to a syntenic location or "*" if cannot be imputed by parsimony

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
    if validRooting(tree, locusMapForRootingD):
        rerootings+=[tree]
    originaltree=frozenset([get_name(left_subtree(tree)),get_name(right_subtree(tree))])
    alltrees.add(originaltree)

    #call the helper function where actually rerooting happens
    helperRerooting(tree, alltrees, rerootings, locusMapForRootingD)
    return rerootings  

def helperRerooting(tree, alltrees, rerootings, locusMapForRootingD):
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
        addTree(tree,newleft,newright, alltrees, rerootings, locusMapForRootingD)
    
        #second we can root on the path to its left
        newleft=(get_name(left_subtree(tree)),left_subtree(left_subtree(tree)),right_subtree(tree))
        newright=right_subtree(left_subtree(tree))

        addTree(tree,newleft,newright, alltrees, rerootings, locusMapForRootingD)
        

    if not is_leaf(right_subtree(tree)):
        #if the right tree has two children
        #first we can root on the path to its right children
        newright=right_subtree(right_subtree(tree))
        newleft=(get_name(right_subtree(tree)),left_subtree(right_subtree(tree)),left_subtree(tree))

        addTree(tree,newleft,newright, alltrees, rerootings, locusMapForRootingD)
    
        #second we can root on the path to its left
        newright=(get_name(right_subtree(tree)),right_subtree(right_subtree(tree)),left_subtree(tree))
        newleft=left_subtree(right_subtree(tree))

        addTree(tree,newleft,newright, alltrees, rerootings, locusMapForRootingD)

def addTree(tree, newleft, newright, alltrees, rerootings, locusMapForRootingD):
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
        if validRooting(newrooting, locusMapForRootingD):
            rerootings+=[newrooting]
        alltrees.add(newtree)
        helperRerooting(newrooting, alltrees, rerootings, locusMapForRootingD)

#### Convert to tree format for DTLOR

def parseTreeForDP(tree, parasite):
    """
    Input: four tuple tree as defined by the trees.py (xenoGI)
    Output: tree format represented by a dictionary as required by the 
            DTLOR_DP.py. Note, if an edge directs to a tip, than the children 
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

#### Another dict format for trees (we should really make a tree object...)

def getTreeDictionary(geneTree,tree_dict):
    """
    Input: tuple tree format 
    Output: dictionary with key as parent node and value as the children tuple
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

def getLeavesInSubtree(startNode, tree_dict):
    """
    returns the list of leaves in the subtree rooted at startNode
    """
    if startNode not in tree_dict: #it is a tip
        return [int(startNode)]    #prepare for addGene which needs integer
    else:
        left, right=tree_dict[startNode]
        leaves=getLeavesInSubtree(left, tree_dict)
        leaves.extend(getLeavesInSubtree(right, tree_dict))
        return leaves
