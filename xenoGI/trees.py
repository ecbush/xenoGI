# functions for loading, manipulating and creating phylogenetic trees
from Bio import Phylo
import glob,os,shutil,subprocess
from multiprocessing import Pool
from .Tree import *
from xenoGI import genomes,fasta

#### Functions working with biopython trees

def prepareTree(bpTree,outGroupTaxaL):
    '''Given a biopython unrooted tree object (e.g. ASTRAL output) root,
name the internal nodes, and strip branch lengths (if
any). outGroupTaxaL provides the outgroups for rooting. Return
resulting tree in four tuple form.

    '''
    # strip confidence values in internal nodes (which mess up our
    # names)
    for n in bpTree.get_nonterminals():
        n.confidence = None
    
    bpTree = rootTree(bpTree,outGroupTaxaL)
    bpTree = nameInternalNodes(bpTree,"s")

    return bpTree

def rootTree(bpTree,outGroupTaxaL):
    '''Given a bio python tree, root it. outGroupTaxaL should consist of
all the species in the desired outgroup clade.'''

    bpTree.root_with_outgroup(*outGroupTaxaL)
    # rooting with outgroups with multiple species seems to leave a
    # trifurcating root.
    if len(bpTree.root)==2:
        return bpTree
    else:
        return makeRootBifurcation(bpTree,outGroupTaxaL)        

def makeRootBifurcation(bpTree,outGroupTaxaL):
    '''Given a bpTree that is not bifurcating at root, if it is
trifurcating convert to bifurcating. Otherwise thow error.'''

    if len(bpTree.root)!=3:
        raise ValueError("This tree has a strange number of branches at the root (not 2 or 3).")

    # find the two in outGroupTaxaL, and the one not in it
    inL=[]
    notIn=None
    for cld in bpTree.root:
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
    outBpTree = Phylo.Newick.Tree(root=c2,rooted=True)

    return outBpTree
    
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

def nameInternalNodes(bpTree,nameStem):
    '''Given a bio python tree without named internal nodes, name
them. This obliterates any names that might have been inside
already (e.g. if confidence values have been put in the name field).'''
    for j,n in enumerate(bpTree.get_nonterminals()):
        n.name=nameStem+str(j)
    return bpTree 

#### Functions for creating gene and species trees

def makeSpeciesTree(paramD,aabrhHardCoreL,genesO):
    '''Create a species tree based on the ortholog sets in aabrhHardCoreL.'''

    ## set up
    # create work dir if it doesn't already exist
    workDir = paramD['makeSpeciesTreeWorkingDir']
    if glob.glob(workDir)==[]:
        os.mkdir(workDir)

    numProcesses = paramD['numProcesses']
    javaPath = paramD['javaPath']
    astralPath = paramD['astralPath']
    astralTreeFN = paramD['astralTreeFN']
    gtFileStem = paramD['aabrhHardCoreGeneTreeFileStem']
    allGtFilePath = os.path.join(workDir,gtFileStem+'*.tre')
    aabrhHardCoreGeneTreesFN = paramD['aabrhHardCoreGeneTreesFN']
    outSpeciesTreeFN = paramD['speciesTreeFN'] # for main output
    outGroupTaxaL = [paramD['outGroup']]
    deleteSpeciesTreeWorkingDir = paramD['deleteSpeciesTreeWorkingDir']
    
    # if tree file already exists, throw error
    if os.path.isfile(outSpeciesTreeFN):
        raise IOError("The tree file " + outSpeciesTreeFN + " already exists.")

    # delete any pre-existing hard core gene trees
    for fn in glob.glob(allGtFilePath):
        os.remove(fn)
    
    ## make gene tree for each aabrh hard Core set
    # add numbering to list
    newAabrhHardCoreL = []
    for orthoNum,orthoT in enumerate(aabrhHardCoreL):
        newAabrhHardCoreL.append((orthoNum,orthoT))
    
    makeGeneTreesFastTree(paramD,True,genesO,workDir,gtFileStem,newAabrhHardCoreL)

    ## run Astral on gene trees

    # concatenate all gene tree files
    with open(aabrhHardCoreGeneTreesFN, 'w') as aabrhHardCoreGeneTreesF:
        for fn in sorted(glob.glob(allGtFilePath)):
            with open(fn) as infile:
                aabrhHardCoreGeneTreesF.write(infile.read())

    # run astral
    try:
        subprocess.check_call([javaPath, '-jar',astralPath, '-i', aabrhHardCoreGeneTreesFN, '-o',astralTreeFN],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    except:
        raise Exception("Astral is giving an error on the file "+aabrhHardCoreGeneTreesFN)
        
    ## load and root the tree.
    speciesRtreeO = Rtree()
    speciesRtreeO.fromNewickFileLoadSpeciesTree(astralTreeFN,outGroupTaxaL)
    with open(outSpeciesTreeFN,"w") as f: # final output
        f.write(speciesRtreeO.toNewickStr()+";\n")

    ## delete the working directory
    if deleteSpeciesTreeWorkingDir:
        shutil.rmtree(workDir)

# GeneRax gene tree making        
def makeGeneTreesGeneRax(paramD,strainHeader,genesO,workDir,gtFileStem,orthoTL):
    '''Given a list of ortho groups, orthoTL, use GeneRax to make a gene
 tree for each and put in workDir. strainHeader is a boolean that if
 True means we put the strain in the header for alignments (so we'll
 end up with that for the names of the tips of the tree). gtFileStem
 gives the stem of the name for the gene tree files. orthoTL is a list
 of (orthoGroupNum,orthoT) where orthoT has the genes in an ortholog
 group.

    '''

    numProcesses = paramD['numProcesses']    
    musclePath = paramD['musclePath']
    geneRaxPath = paramD['geneRaxPath']
    GeneRaxInputListFN = paramD['GeneRaxInputListFN']
    speciesTreeFN = paramD['speciesTreeFN']
    GeneRaxReconcilationModel = paramD['GeneRaxReconcilationModel']
    GeneRaxSearchRadius = paramD['GeneRaxSearchRadius']
    GeneRaxOutputDirN = paramD['GeneRaxOutputDirN']+"-"+gtFileStem
    alignStem = "align"+"-"+gtFileStem
    
    # create alignments
    createAlignments(paramD,alignStem,strainHeader,genesO,workDir,musclePath,numProcesses,orthoTL)

    # create support files for running GeneRax
    createMappingFilesGeneRax(paramD,orthoTL,workDir,genesO)

    # create the file listing alignment and mapping files for each family
    createGeneRaxInputList(paramD,orthoTL,workDir,alignStem)

    # run generax
    subprocess.check_call(["mpiexec", "-np",str(numProcesses),geneRaxPath,"-f",os.path.join(workDir,GeneRaxInputListFN),"-s",speciesTreeFN,"-r",GeneRaxReconcilationModel,"--max-spr-radius",str(GeneRaxSearchRadius),"-p",os.path.join(workDir,GeneRaxOutputDirN)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    
    # copy genetree files out of generax output directory. If generax failed on any, keep track
    failedOrthoSetL = []
    for orthoGroupNum,orthoT in orthoTL:
        grOutputFN = os.path.join(workDir,GeneRaxOutputDirN,"results","family_"+str(orthoGroupNum),"geneTree.newick")
        orthoGroupNumStr = str(orthoGroupNum).zfill(6) # pad w/ 0's
        destGeneTreeFN = os.path.join(workDir,gtFileStem+orthoGroupNumStr+".tre")

        try:
            shutil.copyfile(grOutputFN,destGeneTreeFN)
        except FileNotFoundError:
            failedOrthoSetL.append((orthoGroupNum,orthoT))   

    return failedOrthoSetL
            
            
def createAlignments(paramD,alignStem,strainHeader,genesO,workDir,musclePath,numProcesses,orthoTL):
    '''Create alignments for all the ortho groups in orthoTL using MUSCLE.'''

    # set up argumentL

    argumentL = [([],[],paramD,alignStem,strainHeader,genesO,workDir,musclePath) for i in range(numProcesses)]

    counter = 0
    for orthoGroupNum,orthoT in orthoTL:
        # we use counter for placement in argumentL since at least in
        # the case of single gene families, some will be missing
        orthoGroupNumStr = str(orthoGroupNum).zfill(6) # pad w/ 0's so ls will display in right order
        argumentL[counter%numProcesses][0].append(orthoGroupNumStr)
        argumentL[counter%numProcesses][1].append(orthoT)
        counter += 1
        
    # run
    with Pool(processes=numProcesses) as p:
        for _ in p.imap_unordered(alignOneOrthoTGroup, argumentL):
            pass

    return

def alignOneOrthoTGroup(argT):
    '''Function for aligning one family in a multi processing context.'''
    
    orthoGroupNumStrL,orthoTL,paramD,alignStem,strainHeader,genesO,workDir,musclePath = argT

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

    # align each orthoT in this orthoTL
    for i in range(len(orthoTL)):
        orthoT = orthoTL[i]
        orthoGroupNumStr = orthoGroupNumStrL[i]

        # get temp and output align file names
        inTempProtFN=os.path.join(workDir,"tempProt-"+orthoGroupNumStr+".fa")
        outAlignFN=os.path.join(workDir,alignStem+"-"+orthoGroupNumStr+".afa")
        
        if len(orthoT) <2: # if there's only one seq, make that the align file
            writeFasta(outAlignFN,orthoT,strainHeader,genesO,protSeqD)
        else:
            alignOneOrthoT(orthoT,strainHeader,musclePath,inTempProtFN,outAlignFN,protSeqD,dnaSeqD,genesO)
            os.remove(inTempProtFN) # remove unaligned prot file

def alignOneOrthoT(orthoT,strainHeader,musclePath,inProtFN,outAlignFN,protSeqD,dnaSeqD,genesO):
    '''Given genes in a single ortholog set, align them with muscle. If
dnaSeqD is empty, uses protein only.'''

    # write prots we want to align to temp file
    writeFasta(inProtFN,orthoT,strainHeader,genesO,protSeqD)

    # align proteins (muscle 5 syntax here. only 1 thread each since we're doing many)
    retCode = subprocess.call([musclePath,'-threads', '1', '-align' ,inProtFN, '-output', outAlignFN],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

    if retCode != 0:
        raise Exception("Alignment failed for "+inProtFN)
    
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
        
def createMappingFilesGeneRax(paramD,orthoTL,workDir,genesO):
    '''Create the set of mapping files that shows which species each gene is in.'''
    
    mappingFileStem = paramD['GeneRaxMappingFileStem']
    mappingFileExt = paramD['GeneRaxMappingFileExt']

    for orthoGroupNum,orthoT in orthoTL:
        orthoGroupNumStr = str(orthoGroupNum).zfill(6) # pad w/ 0's so ls will display in right order
        mappingFN = os.path.join(workDir,mappingFileStem + orthoGroupNumStr + mappingFileExt)
        
        with open(mappingFN,'w') as mappingF:
            for geneNum in orthoT:
                strainName = genesO.numToStrainName(geneNum)
                mappingF.write(str(geneNum)+" "+strainName+"\n")

def createGeneRaxInputList(paramD,orthoTL,workDir,alignStem):
    '''Create the geneRaxInputList.txt file which asigns alignment and
mapping files to each family.'''

    mappingFileStem = paramD['GeneRaxMappingFileStem']
    mappingFileExt = paramD['GeneRaxMappingFileExt']
    GeneRaxInputListFN = paramD['GeneRaxInputListFN']
    
    if paramD['dnaBasedGeneTrees']:
        subModel = paramD['GeneRaxDNASubstModel']
    else:
        subModel = paramD['GeneRaxProtSubstModel']
        
    fn = os.path.join(workDir,GeneRaxInputListFN)
    
    with open(fn,'w') as f:

        f.write("[FAMILIES] # this is a comment\n")

        for orthoGroupNum,orthoT in orthoTL:
            orthoGroupNumStr = str(orthoGroupNum).zfill(6) # pad w/ 0's so ls will display in right order
            # we will run GeneRax from workDir, so don't need to include in paths
            mappingFN = os.path.join(workDir,mappingFileStem + orthoGroupNumStr + mappingFileExt)
            alignFN=os.path.join(workDir,alignStem+"-"+orthoGroupNumStr+".afa")

            f.write("- family_"+str(orthoGroupNum)+"\n")
            f.write("alignment = "+alignFN+"\n")
            f.write("mapping = "+mappingFN+"\n")
            f.write("subst_model = "+subModel+"\n")

# FastTree gene tree making        
def makeGeneTreesFastTree(paramD,strainHeader,genesO,workDir,gtFileStem,orthoTL):
    '''Given a list of ortho groups, orthoTL, use FastTree to make a gene
 tree for each and put in workDir. strainHeader is a boolean that if
 True means we put the strain in the header for alignments (so we'll
 end up with that for the names of the tips of the tree). gtFileStem
 gives the stem of the name for the gene tree files. orthoTL is a list
 of (orthoGroupNum,orthoT) where orthoT has the genes in an ortholog
 group.

    '''
    alignStem = "align"+"-"+gtFileStem
    numProcesses = paramD['numProcesses']    
    musclePath = paramD['musclePath']
    fastTreePath = paramD['fastTreePath']

    ## make alignments
    createAlignments(paramD,alignStem,strainHeader,genesO,workDir,musclePath,numProcesses,orthoTL)
    
    ## make gene trees
    argumentL = [([],paramD,workDir,gtFileStem,fastTreePath) for i in range(numProcesses)]

    # set up argumentL
    counter = 0
    for orthoGroupNum,orthoT in orthoTL:
        # we use counter for placememnt in argumentL since at least in
        # the case of single gene families, some will be missing
        orthoGroupNumStr = str(orthoGroupNum).zfill(6) # pad w/ 0's so ls will display in right order
        argumentL[counter%numProcesses][0].append(orthoGroupNumStr)
        counter += 1
        
    # run
    with Pool(processes=numProcesses) as p:
        for _ in p.imap_unordered(makeOneGeneTreeGroupFastTree, argumentL):
            pass

    return
        
def makeOneGeneTreeGroupFastTree(argT):
    '''Wrapper for multiprocessing. Given a group of orthoT's to work
    on, calls makeOneGeneTreeFastTree on each.'''
    
    orthoGroupNumStrL,paramD,workDir,gtFileStem,fastTreePath = argT

    dnaBasedGeneTrees = paramD['dnaBasedGeneTrees']

    # make trees
    for orthoGroupNumStr in orthoGroupNumStrL:
        makeOneGeneTreeFastTree(gtFileStem,orthoGroupNumStr,workDir,dnaBasedGeneTrees,fastTreePath)
    
def makeOneGeneTreeFastTree(gtFileStem,orthoGroupNumStr,workDir,dnaBasedGeneTrees,fastTreePath):
    ''' Makes one gene tree from an ortho list. If dnaBasedGeneTrees is False, uses protein only.'''
 
    # get align file names
    alignStem = "align"+"-"+gtFileStem
    outAlignFN=os.path.join(workDir,alignStem+"-"+orthoGroupNumStr+".afa")

    ## make gene tree
    geneTreeFN = os.path.join(workDir,gtFileStem+orthoGroupNumStr+".tre")
    if dnaBasedGeneTrees == False:
        # using protein
        subprocess.check_call([fastTreePath, '-out',geneTreeFN,outAlignFN],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    else:
        # using dna
        subprocess.check_call([fastTreePath,'-gtr','-nt','-out',geneTreeFN,outAlignFN],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    
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

