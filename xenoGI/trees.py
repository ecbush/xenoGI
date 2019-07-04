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

    numThreads = paramD['numThreads']    
    musclePath = paramD['musclePath']
    fastTreePath = paramD['fastTreePath']
    javaPath = paramD['javaPath']
    astralPath = paramD['astralPath']
    astralTreeFN = paramD['astralTreeFN']
    gtFileStem = 'aabrhHardCoreFam'
    allGtFilePath = os.path.join(workDir,gtFileStem+'*.tre')
    aabrhHardCoreGeneTreesFN = paramD['aabrhHardCoreGeneTreesFN']
    outSpeciesTreeFN = paramD['treeFN'] # for main output
    deleteSpeciesTreeWorkingDir = paramD['deleteSpeciesTreeWorkingDir']
    outGroupTaxaL = paramD['outGroupTaxaL']
    
    # if tree file already exists, throw error
    if os.path.isfile(paramD['treeFN']):
        raise IOError("The tree file " + paramD['treeFN'] + " already exists.")
    
    # load protein and dna sequences
    protSeqD=genomes.loadSeq(paramD, '_prot.fa')

    if paramD['dnaBasedSpeciesTree'] == True:
        dnaSeqD = genomes.loadSeq(paramD, '_dna.fa')
    else:
        dnaSeqD = {}

    ## make gene trees
    argumentL = [([],[],genesO,protSeqD,dnaSeqD,workDir,gtFileStem,musclePath,fastTreePath) for i in range(numThreads)]

    # set up argumentL
    orthoGroupNum = 0
    for orthoT in aabrhHardCoreL:
        orthoGroupNumStr = str(orthoGroupNum).zfill(6) # pad with 0's so ls will display in right order
        argumentL[orthoGroupNum%numThreads][0].append(orthoGroupNumStr)
        argumentL[orthoGroupNum%numThreads][1].append(orthoT)
        orthoGroupNum+=1

    # run
    with Pool(processes=numThreads) as p:
        for _ in p.imap_unordered(makeOneGeneTreeGroup, argumentL):
            pass

    ## run Astral on gene trees

    # concatenate all gene tree files
    with open(aabrhHardCoreGeneTreesFN, 'w') as aabrhHardCoreGeneTreesF:
        for fn in sorted(glob.glob(allGtFilePath)):
            with open(fn) as infile:
                aabrhHardCoreGeneTreesF.write(infile.read())
    
    subprocess.call([javaPath, '-jar',astralPath, '-i', aabrhHardCoreGeneTreesFN, '-o',astralTreeFN],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    ## load and root the tree.
    bpTree = Phylo.read(astralTreeFN, 'newick')
    tupleTree=prepareTree(bpTree,outGroupTaxaL)
    writeTreeNoBrLen(tupleTree,outSpeciesTreeFN) # final output to 
    
    ## delete the working directory
    if deleteSpeciesTreeWorkingDir:
        shutil.rmtree(workDir)

def makeOneGeneTreeGroup(argT):
    '''Wrapper for multiprocessing. Given a group of orthoT's to work on,
calls makeOneGeneTree on each.'''
    
    orthoGroupNumStrL,orthoTL,genesO,protSeqD,dnaSeqD,workDir,gtFileStem,musclePath,fastTreePath = argT

    for i in range(len(orthoTL)):
        orthoGroupNumStr = orthoGroupNumStrL[i]
        orthoT = orthoTL[i]
        makeOneGeneTree(orthoGroupNumStr,orthoT,genesO,protSeqD,dnaSeqD,workDir,gtFileStem,musclePath,fastTreePath)
    
def makeOneGeneTree(orthoGroupNumStr,orthoT,genesO,protSeqD,dnaSeqD,workDir,gtFileStem,musclePath,fastTreePath):
    ''' Makes one gene tree from an ortho list. If dnaSeqD is empty, uses protein only.'''
 
    # get temp align file names
    inTempProtFN=os.path.join(workDir,"tempProt"+orthoGroupNumStr+".fa")
    outAlignFN=os.path.join(workDir,"align"+orthoGroupNumStr+".afa")

    # align
    alignOneOrthoT(orthoT,musclePath,inTempProtFN,outAlignFN,protSeqD,dnaSeqD,genesO)
    os.remove(inTempProtFN) # remove unaligned prot file
    
    # make gene tree
    geneTreeFN = os.path.join(workDir,gtFileStem+orthoGroupNumStr+".tre")
    if dnaSeqD == {}:
        # using protein
        subprocess.call([fastTreePath, '-out',geneTreeFN,outAlignFN],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    else:
        # using dna
        subprocess.call([fastTreePath,'-gtr','-nt','-out',geneTreeFN,outAlignFN],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

def alignOneOrthoT(orthoT,musclePath,inProtFN,outAlignFN,protSeqD,dnaSeqD,genesO):
    '''Given genes in a single ortholog set, align them with muscle. If
dnaSeqD is empty, uses protein only.'''

    # write prots we want to align to temp file
    with open(inProtFN,"w") as tempf:
        writeSeqBlock(tempf,orthoT,genesO,protSeqD)

    # align proteins
    subprocess.call([musclePath, '-in' ,inProtFN, '-out', outAlignFN],stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    if dnaSeqD != {}:
        # back align to get dna alignment, overwriting protein alignment.
        protAlignL = []
        for hd,alignedProtSeq in fasta.load(outAlignFN):
            protAlignL.append((int(hd.rstrip().split()[1]),alignedProtSeq))
        backAlign(outAlignFN,protAlignL,dnaSeqD,genesO)

    return
    
def writeSeqBlock(f,orthoT,genesO,seqD):
    '''writes a multifasta block. seqs are specified in
orthoT, and obtained from seqD.'''
    for geneNum in orthoT:
        strainName = genesO.numToStrainName(geneNum)
        f.write(">" + strainName + ' '  + str(geneNum) + "\n")
        f.write(seqD[geneNum] + "\n")
        f.write("\n")

def backAlign(outAlignFN,protAlignL,dnaSeqD,genesO):
    '''Prints the nucleotide alignments given by one block of protein
    alignments.'''
    with open(outAlignFN,'w') as outf:
        printBlock=''
        for geneNum,protSeq in protAlignL:
            dnaSeq = dnaSeqD[geneNum]
            lenProtein = len(protSeq) - protSeq.count('-')
            #protLength(protSeq)
            # if the protein and sequence don't match in length, throw it out 
            # we have to check two cases because some proteins include the stop
            # codon and some don't
            if lenProtein * 3 != len(dnaSeq) and (lenProtein + 1) * 3 != len(dnaSeq):
                raise IndexError("Lengths of dna and protein do not correspond.")

            strainName = genesO.numToStrainName(geneNum)
            # header
            printBlock += ">" + strainName + ' ' + str(geneNum) + "\n"
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
