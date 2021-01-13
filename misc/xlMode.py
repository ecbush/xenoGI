import sys,random,os,glob,numpy,ctypes,math
sys.path.insert(0,os.path.join(sys.path[0],'../xenoGI'))
from xenoGI import xenoGI,blast,trees,scores,genomes,Score,families,islands,analysis
from Bio import Phylo

random.seed(42)

###### Obtain core ortho sets wrapper ######
        
def obtainCoreOrthoSetsWrapper(paramD):
    """Wrapper to obtain an initial set of orthologs that can be used to
make a species tree. We 1. Take a random sample of genomes and get
sets of all around best reciprocal hit orthologs for these. 2. Blast
the aabrh gene set vs. all the genomes. 3. Identify those sets of
aabrh orthologs that have hits in all the genomes and write these to
file."""
    ## pick a random subset of genomes and get aabrh sets for these
    allGenomesFastaL=glob.glob(paramD['fastaFilePath'])
    allGenomesFastaL = [fastaFN.rstrip() for fastaFN in allGenomesFastaL if '_prot' in fastaFN]
    
    # if sample file already exists, read it in
    if os.path.isfile(paramD['genomesInRandomSampleFN']):
        sampleGenomesFastaL = [fn.rstrip() for fn in open(paramD['genomesInRandomSampleFN'])]
    else:
        # take random sample and write to file
        sampleGenomesFastaL = random.sample(allGenomesFastaL,paramD['numGenomesInRandomSample']-1)
        # add outgroup
        for fn in allGenomesFastaL:
            if paramD['outGroup'] in fn:
                sampleGenomesFastaL.append(fn)
                break
        # save sample genomes to file
        f = open(paramD['genomesInRandomSampleFN'], 'w')
        for fn in sampleGenomesFastaL:
            f.write(fn+"\n")
        f.close()
        
    # run blast to get blast files needed to create the aabrhL
    blast.runBlast(sampleGenomesFastaL,sampleGenomesFastaL,paramD)
    sampleGenomesStrainNamesL = []
    for genomeFastaFN in sampleGenomesFastaL:
        strainName = os.path.split(genomeFastaFN)[-1] # strain + extension
        strainName = os.path.splitext(strainName)[0] # strain only
        strainName = strainName.split('_prot')[0]
        sampleGenomesStrainNamesL.append(strainName)
    allGenomesStrainNamesL = []
    for genomeFastaFN in allGenomesFastaL:
        strainName = os.path.split(genomeFastaFN)[-1]
        strainName = os.path.splitext(strainName)[0]
        strainName = strainName.split('_prot')[0]
        allGenomesStrainNamesL.append(strainName)

    # when testing using existing fasta files. uncomment out for new genomes/empty fasta file
    if 'randomSampleAabrh' in allGenomesStrainNamesL:
        allGenomesStrainNamesL.remove('randomSampleAabrh')

    # get list of all around best reciprocal hits for this sample
    randomSampleAabrhL = scores.createAabrhL(paramD['blastFilePath'],sampleGenomesStrainNamesL,paramD['evalueThresh'],paramD['alignCoverThresh'],paramD['percIdentThresh'],paramD['randomSampleAabrhFN'])
    
    # write all to a single file
    fastaDir = os.path.split(paramD['fastaFilePath'])[0]
    randomSampleAabrhFastaFilePath = os.path.join(fastaDir,paramD['randomSampleAabrhFastaFN'])
    createAabrhFasta(sampleGenomesFastaL,randomSampleAabrhFastaFilePath,randomSampleAabrhL,paramD)

    # now we blast randomSampleAabrhFasta vs. all genomes
    blast.runBlast([randomSampleAabrhFastaFilePath],allGenomesFastaL, paramD)

    # create list of lists of conserved core gene set
    orthoL = [[] for aabrhT in randomSampleAabrhL]
    blastDir = os.path.split(paramD['blastFilePath'])[0]
    orthoL = getCoreGenes(allGenomesStrainNamesL,blastDir,randomSampleAabrhL, orthoL, paramD)	

    # get rid of Nones in list
    orthoL = [i for i in orthoL if i]

    # write into file
    f = open(paramD['allStrainCoreOrthosFN'],"w")
    for ortho in orthoL:
        for gene in ortho:
            f.write(str(gene) + '\t')
        f.write('\n')
    f.close()
    return orthoL

def createAabrhFasta(sampleGenomesFastaL,randomSampleAabrhFastaFN,aabrhL, paramD):
    '''If it doesn't already exist, create a fasta file containing all the
proteins in the initial random sample aabrh set. If
randomSampleAabrhFastaFN does exist, we simply leave it.
    '''
    if not os.path.isfile(paramD['randomSampleAabrhFastaFN']):
        # load sequence dict, but only for genes in aabrhL
        randSampGenesS = set()
        for orthoT in aabrhL:
            randSampGenesS.update(orthoT)
        seqD = genomes.loadSeq(paramD, '_prot.fa',randSampGenesS)
        # write fasta
        f=open(randomSampleAabrhFastaFN,'w')
        for orthoT in aabrhL:
            for gene in orthoT:
                writeFastaLine(f,gene,seqD)
        f.close()
    return

def writeFastaLine(f,gene,seqD):
    '''Write a fasta line of gene to file handle f.'''
    f.write(">"+str(gene)+"_rsamp"+"\n")
    f.write(seqD[gene]+"\n")

def getCoreGenes(allGenomesStrainNamesL,blastDir,randomSampleAabrhL,orthoL,paramD):
    # for every strain in all genomes, we find all blast hits with the sampled group. 
    for strain in allGenomesStrainNamesL:
        randomSampleAabrhFastaStem = paramD['randomSampleAabrhFastaFN'].split(".fa")[0]
        blastFN = os.path.join(blastDir,randomSampleAabrhFastaStem+'_-VS-_'+strain+'.out')
        strainHitsD = scores.getHits(blastFN,paramD['evalueThresh'],paramD['alignCoverThresh'], paramD['percIdentThresh'])
        # now for every set, we check if the best hit of each gene is the same gene 
        # if so, we add that gene to the orthoL and move to the next genome. else, remove the set
        for aabrhInd in range(len(randomSampleAabrhL)):
            # check if we've already gotten rid of the set
            if orthoL[aabrhInd] != None:
                gene = randomSampleAabrhL[aabrhInd][0]
                # check if gene is in dictionary, if not get rid of set. 
                if gene in strainHitsD:
                    strainGene = strainHitsD[gene]
                    # check if all genes have same hit, if so add gene to the set
                    if orthoChecker(randomSampleAabrhL[aabrhInd], strainHitsD, strainGene):
                        orthoL[aabrhInd].append(strainGene)
                    else:
                        orthoL[aabrhInd] = None
                else:
                    orthoL[aabrhInd] = None
            else:
                continue
    return orthoL

def orthoChecker(ortho, strainHitsD, strainGene):
    """ sub-function for obtainCoreGenesWrapper function to check whether all the genes in 
the tuple have the same hit in a genome"""
    # check if gene is in dictionary, if so check if every gene gives the same hit. if so, return true
    for gene in ortho:
        if gene not in strainHitsD:
            return False
        else:
            if strainHitsD[gene] == strainGene:
               return True
            else:
               return False


###### Trim tree ######

def trimTree(paramD):
    '''Load a tree and trim it to trimLeafNum using the treemer algorithm.'''

    # parameters
    trimLeafNum = paramD['trimLeafNum']
    astralTreeFN = paramD['astralTreeFN']
    outGroupL = [paramD['outGroup']]
    userSpecifiedStrainsFileName = paramD['userSpecifiedStrainsFileName']
    userSpecifiedStrainsL = []
    if userSpecifiedStrainsFileName != None:
        with open(userSpecifiedStrainsFileName,'r') as f:
            while True:
                s=f.readline()
                if s=='':
                    break
                userSpecifiedStrainsL.append(s.rstrip())
    protectedList = list(set(userSpecifiedStrainsL + outGroupL))
    scaffoldTreeFN = paramD['scaffoldTreeFN']
    
    # load tree
    treeToTrimO = trees.Rtree()
    treeToTrimO.fromNewickFileLoadSpeciesTree(astralTreeFN,outGroupL,True)
    prepareBL(treeToTrimO)

    # trim
    startNumLeaves = treeToTrimO.leafCount()
    numIters = startNumLeaves - trimLeafNum
    for i in range(numIters):
        redundList = []
        pairDict = {}
        makeDictForTreeO(treeToTrimO, pairDict, redundList)
        leaf = pickToPruneForTreeO(pairDict, protectedList)
        treeToTrimO = pruneLeafForTreeO(leaf, treeToTrimO)

    # write to file
    with open(scaffoldTreeFN,"w") as f: # final output
        f.write(treeToTrimO.toNewickStr()+";\n")
    

##### Recreating Treemmer for TreeO #####

def getDistForTreeO(node1, node2, treeO):
    ''' given 2 nodes, returns the branchLen between them '''
    branchPair = (node1, node2)
    if branchPair in treeO.branchLenD:
        dist = treeO.branchLenD[branchPair]
    else:
        dist = treeO.branchLenD[(branchPair[1], branchPair[0])]
    return dist

def makeDictForTreeO(treeO, pairDict, redundList):
    ''' calls getNeighborForTreeO on all of the leaves '''
    for leaf in treeO.leaves():
        parentNode = treeO.getParent(leaf)
        if leaf not in redundList and parentNode != treeO.rootNode:
            getNeighborForTreeO(treeO, pairDict,redundList, leaf, parentNode)

def getNeighborForTreeO(treeO, pairDict, redundList, leaf, parentNode):
    ''' finds nearest neighbor of each leaf '''
    leafDist = getDistForTreeO(leaf, parentNode, treeO)
    leftChild, rightChild = treeO.children(parentNode)
    if treeO.isLeaf(leftChild) and leftChild != leaf:
        leftDist = getDistForTreeO(leftChild, parentNode, treeO)
        dist = leafDist + leftDist
        pairDict[(leaf,leftChild)] = dist
        redundList.append(leftChild)
    elif treeO.isLeaf(rightChild) and rightChild != leaf:
        rightDist = getDistForTreeO(rightChild, parentNode, treeO)
        dist = leafDist + rightDist
        pairDict[(leaf,rightChild)] = dist
        redundList.append(rightChild)
    else:
        sisDict = getSisterForTreeO(leaf, parentNode, treeO, leafDist)
        multiDict = getMultiForTreeO(leaf, parentNode, treeO, leafDist)
        pairDict.update(sisDict)
        pairDict.update(multiDict)

def getSisterForTreeO(leaf, parentNode, treeO, leafDist):
    ''' gets the sister of a leaf in a TreeO '''
    tempDict = {}
    leftChild, rightChild = treeO.children(parentNode)
    if not treeO.isLeaf(leftChild):
        leftDist = getDistForTreeO(leftChild, parentNode, treeO)
        leftLeft, leftRight = treeO.children(leftChild)
        if treeO.isLeaf(leftLeft) and not treeO.isLeaf(leftRight):
            leftLeftDist = getDistForTreeO(leftLeft, leftChild, treeO)
            dist = leafDist + leftDist + leftLeftDist
            tempDict[(leaf,leftLeft)] = dist
        if treeO.isLeaf(leftRight) and not treeO.isLeaf(leftLeft):
            leftRightDist = getDistForTreeO(leftRight, leftChild, treeO)
            dist = leafDist + leftDist + leftRightDist
            tempDict[(leaf,leftRight)] = dist
    else:
        rightDist = getDistForTreeO(rightChild, parentNode, treeO)
        rightLeft, rightRight = treeO.children(rightChild)
        if treeO.isLeaf(rightLeft) and not treeO.isLeaf(rightRight):
            rightLeftDist = getDistForTreeO(rightLeft, rightChild, treeO)
            dist = leafDist + rightDist + rightLeftDist
            tempDict[(leaf,rightLeft)] = dist
        if treeO.isLeaf(rightRight) and not treeO.isLeaf(rightLeft):
            rightRightDist = getDistForTreeO(rightRight, rightChild, treeO)
            dist = leafDist + rightDist + rightRightDist
            tempDict[(leaf,rightRight)] = dist
    return tempDict

def getMultiForTreeO(leaf, parentNode, treeO, leafDist):
    ''' gets the multi a leaf in a TreeO '''
    tempDict = {}
    grandparentNode = treeO.getParent(parentNode)
    if grandparentNode != None:
        parentDist = getDistForTreeO(parentNode, grandparentNode, treeO)
        if grandparentNode != treeO.rootNode:
            leftPar, rightPar = treeO.children(grandparentNode)
            if treeO.isLeaf(leftPar):
                leftDist = getDistForTreeO(leftPar, grandparentNode, treeO)
                dist = leafDist + parentDist + leftDist
                tempDict[(leaf,leftPar)] = dist
            if treeO.isLeaf(rightPar):
                rightDist = getDistForTreeO(rightPar, grandparentNode, treeO)
                dist = leafDist + parentDist + rightDist
                tempDict[(leaf,rightPar)] = dist
    return tempDict

def pickToPruneForTreeO(pairDict,protectedList):
    ''' goes through dictionary and picks which leaf should be pruned '''
    minPair = (min(pairDict, key = pairDict.get))
    if minPair[0] in protectedList and minPair[1] in protectedList:
        pairDict.pop(minPair)
        pruneLeaf = pickToPruneForTreeO(pairDict, protectedList)
        return pruneLeaf
    elif minPair[0] in protectedList:
        pruneLeaf = minPair[1]
    elif minPair[1] in protectedList:
        pruneLeaf = minPair[0]
    else:
        pruneLeaf = random.choice(minPair)
    return pruneLeaf

def pruneLeafForTreeO(leaf, treeO):
    ''' prunes the leaf in the tree '''
    parentNode = treeO.getParent(leaf)
    grandparentNode = treeO.getParent(parentNode)
    leftChild, rightChild = treeO.children(parentNode)
    if leftChild == leaf:
        newNodeConnectD, newBranchLenD = fixDictionaries(treeO, leaf, rightChild, parentNode, grandparentNode)
    else:
        newNodeConnectD, newBranchLenD = fixDictionaries(treeO, leaf, leftChild, parentNode, grandparentNode)
    trimmedTreeO = trees.Rtree()
    trimmedTreeO.populateAttributes(newNodeConnectD,treeO.rootNode,newBranchLenD)
    return trimmedTreeO

def fixDictionaries(treeO, remover, keeper, parentNode, grandparentNode):
    newNodeConnectD = treeO.nodeConnectD
    newBranchLenD = treeO.branchLenD
    # update nodeConnectD
    # connect grandparent to keeper
    currentGPConnects = newNodeConnectD[grandparentNode]
    newGPConnects = ()
    for connect in currentGPConnects:
        if connect != parentNode:
            newGPConnects += (connect,)
    newGPConnects += (keeper,)
    newNodeConnectD[grandparentNode] = newGPConnects
    # connect keeper to grandparent
    keeperTup = tuple([connected for connected in treeO.nodeConnectD[keeper] if connected != parentNode])
    newNodeConnectD[keeper] =  (grandparentNode,) + keeperTup 
    # remove parent and remover from dict
    newNodeConnectD.pop(parentNode)
    newNodeConnectD.pop(remover)

    # update branchLenD
    parentDist = getDistForTreeO(parentNode, grandparentNode, treeO)
    keepDist = getDistForTreeO(parentNode, keeper, treeO)
    dist = parentDist + keepDist
    # remove grandparent, parent branch
    gpPbranchPair = (grandparentNode, parentNode)
    if gpPbranchPair in newBranchLenD:
        newBranchLenD.pop(gpPbranchPair)
    else:
        newBranchLenD.pop((gpPbranchPair[1], gpPbranchPair[0]))
    # remove parent, remover branch
    pRbranchPair = (parentNode, remover)
    if pRbranchPair in newBranchLenD:
        newBranchLenD.pop(pRbranchPair)
    else:
        newBranchLenD.pop((pRbranchPair[1], pRbranchPair[0]))
    # add grandparent, keeper branch
    newBranchLenD[(grandparentNode, keeper)] = dist

    return newNodeConnectD, newBranchLenD
        
def prepareBL(treeO):
    '''Adjustment to handle astral trees.'''
    for key, value in treeO.branchLenD.items():
        if treeO.branchLenD[key] == None:
            treeO.branchLenD[key] = 1 # fix astral None branch lengths

def prune(allStrainsTreeO, scafStrainsL, scaffoldTreeFN):
    treeToTrimO = allStrainsTreeO
    startNumLeaves = treeToTrimO.leafCount()
    numIters = startNumLeaves - len(scafStrainsL)
    for i in range(numIters):
        redundList = []
        pairDict = {}
        makeDictForTreeO(treeToTrimO, pairDict, redundList)
        leaf = pickToPruneForTreeO(pairDict, scafStrainsL)
        treeToTrimO = pruneLeafForTreeO(leaf, treeToTrimO)
    
    # write to file
    with open(scaffoldTreeFN,"w") as f: # final output
        f.write(treeToTrimO.toNewickStr()+";\n")
    return treeToTrimO

def pruneGeneTree(geneTreeO, numGenes):
    ''' method to trim geneTrees to get repGenes '''
    treeToTrimO = geneTreeO
    startNumLeaves = treeToTrimO.leafCount()
    numIters = startNumLeaves - numGenes
    for i in range(numIters):
        redundList = []
        pairDict = {}
        makeDictForTreeO(treeToTrimO, pairDict, redundList)
        leaf = pickToPruneForTreeO(pairDict,[])
        treeToTrimO = pruneLeafForTreeO(leaf, treeToTrimO)
    return list(treeToTrimO.leaves())


###### Map all genes to scaffold, improve scaffold ######

def scaffoldWrapper(paramD):
    '''Wrapper function that runs xenoGI on the initial scaffold tree,
maps all genes onto this, chooses additional strains to add to the
scaffold, and finally maps all genes again.'''

    ## parameters
    allStrainsTreeFN = paramD['astralTreeFN']
    scaffoldTreeFN = paramD['scaffoldTreeFN']
    fastaDir = os.path.split(paramD['fastaFilePath'])[0]
    blastDir = os.path.split(paramD['blastFilePath'])[0]
    scaffoldFamilyRepGenesFastaFN = os.path.join(fastaDir,paramD['scaffoldFamilyRepGenesFastaFN'])
    scaffoldFamilyRepGenesNumPerFamily = paramD['scaffoldFamilyRepGenesNumPerFamily']
    strainInfoFN = paramD['strainInfoFN']
    
    ## Load a few useful data structues
    genesO = genomes.genes(paramD['geneInfoFN'])
    outGroupL = [paramD['outGroup']]
    scaffoldTree = trees.Rtree()
    scaffoldTree.fromNewickFileLoadSpeciesTree(scaffoldTreeFN)
    allStrainsT = xenoGI.readStrainInfoFN(strainInfoFN)


    ## Add in the user specified strains (Consider adding functionality to specify strains without choosing all)
    '''
    scaffoldStrainsL = list(scaffoldTree.leaves())
    userSpecifiedStrainsFileName = paramD['userSpecifiedStrainsFileName']
    userSpecifiedStrainsL = []
    if userSpecifiedStrainsFileName != None:
        f=open(userSpecifiedStrainsFileName,'r')
        while True:
            s=f.readline()
            if s=='':
                break
            userSpecifiedStrainsL.append(s)
    
        
    newScaffoldStrainsL = list(set(scaffoldStrainsL) | set(userSpecifiedStrainsL))
    allStrainsTree = trees.Rtree()
    allStrainsTree.fromNewickFileLoadSpeciesTree(allStrainsTreeFN,outGroupL,True)
    prepareBL(allStrainsTree)
    scaffoldTree = prune(allStrainsTree,newScaffoldStrainsL, scaffoldTreeFN)
    '''
    
    ## run xenoGI
    print("  xenoGI first run",file=sys.stderr)
    initialFamiliesO,originFamiliesO = makeFamiliesScaffold(paramD,scaffoldTree,fastaDir,blastDir,genesO)
    
    ## blast representative genes from each family vs. all strains
    print("  blast rep vs. all",file=sys.stderr)
    blastRepGenesVsScaffold(paramD,originFamiliesO,scaffoldTree,scaffoldFamilyRepGenesFastaFN,scaffoldFamilyRepGenesNumPerFamily,genesO,allStrainsT,fastaDir)

    ## map genes using blast output
    print("  map all genes",file=sys.stderr)
    locFamNumAr,unmappedVal = mapAllStrainsGenes(paramD,genesO.numGenes,originFamiliesO,allStrainsT,blastDir)
    mappedCt,unmappedCt,mappedProp,unmappedProp = arrayCounts(locFamNumAr,unmappedVal)

    rowsL = [['mappedCt','unmappedCt','mappedProp','unmappedProp']]
    rowsL.append([str(mappedCt),str(unmappedCt),format(mappedProp,".3f"),format(unmappedProp,".3f")])
    analysis.printTable(rowsL,indent=2,fileF=sys.stderr)

    # temp
    #locFamNumAr,unmappedVal = readAllGenesMapToScaffoldAr(paramD['allGenesMapToScaffoldFN'])

    ## get unmapped genes and put them in a db, and blast it vs. itself
    print("  blast unmapped",file=sys.stderr)
    unMappedGenesFileL = writeUnmappedToFile(paramD,locFamNumAr,unmappedVal,fastaDir)
    blast.runBlast(unMappedGenesFileL,unMappedGenesFileL,paramD)

    ## pick additional strains and add to scaffold tree
    print("  pick additional strains",file=sys.stderr)
    scaffoldStrainsL = list(scaffoldTree.leaves())
    strainsToAddL = pickAdditionalStrains(allStrainsT,scaffoldStrainsL,blastDir,unMappedGenesFileL,paramD['evalueThresh'],paramD['xlMapAlignCoverThresh'],genesO,paramD['numStrainsToAddToScaffold'],paramD['percIdentThresh'])
    allStrainsTree = trees.Rtree()
    allStrainsTree.fromNewickFileLoadSpeciesTree(allStrainsTreeFN,outGroupL,True)
    prepareBL(allStrainsTree)
    scaffoldTree = prune(allStrainsTree,scaffoldStrainsL+strainsToAddL, scaffoldTreeFN) # new scaffold
    scaffoldStrainsL = list(scaffoldTree.leaves())
    
    ## re-run xenoGI on larger scaffold
    print("  xenoGI second run",file=sys.stderr)
    rootFocalClade = "".join([x for x in scaffoldTree.children(scaffoldTree.rootNode) if x != outGroupL[0]]) # get child that is not outgroup
    geneOrderD=genomes.createGeneOrderD(paramD['geneOrderFN'],scaffoldStrainsL)
    subtreeD=scaffoldTree.createSubtreeD()
    familiesO = makeFamiliesScaffold(paramD,scaffoldTree,fastaDir,blastDir,genesO)
    originFamiliesO = familiesO[1]
    # also do islands this time
    with open(paramD['islandFormationSummaryFN'],'w') as islandFormationSummaryF:
        locIslByNodeD = islands.makeLocusIslands(geneOrderD,subtreeD,scaffoldTree,paramD,originFamiliesO,rootFocalClade,islandFormationSummaryF)

    ## blast representative genes from each family vs. all strains
    print("  blast rep vs. all again",file=sys.stderr)
    blastRepGenesVsScaffold(paramD,originFamiliesO,scaffoldTree,scaffoldFamilyRepGenesFastaFN,scaffoldFamilyRepGenesNumPerFamily,genesO,allStrainsT,fastaDir)
        
    ## final mapping
    print("  map all genes again",file=sys.stderr)
    locFamNumAr,unmappedVal = mapAllStrainsGenes(paramD,genesO.numGenes,originFamiliesO,allStrainsT,blastDir)
    mappedCt,unmappedCt,mappedProp,unmappedProp = arrayCounts(locFamNumAr,unmappedVal)
    rowsL = [['mappedCt','unmappedCt','mappedProp','unmappedProp']]
    rowsL.append([str(mappedCt),str(unmappedCt),format(mappedProp,".3f"),format(unmappedProp,".3f")])
    analysis.printTable(rowsL,indent=2,fileF=sys.stderr)
    
    return

def makeFamiliesScaffold(paramD,scaffoldTree,fastaDir,blastDir,genesO):
    '''Run regular xenoGI on a scaffold tree.'''

    # some useful stuff
    scaffoldStrainL = list(scaffoldTree.leaves())
    geneOrderD = genomes.createGeneOrderD(paramD['geneOrderFN'], scaffoldStrainL)
    subtreeD = scaffoldTree.createSubtreeD
    
    # blast scaffold strains
    fastaToBlastL = []
    for strain in scaffoldStrainL:
        fastaToBlastL.append(os.path.join(fastaDir, strain+'_prot.fa'))
        
    blast.runBlast(fastaToBlastL,fastaToBlastL,paramD)

    blastFnL = []
    for strain1 in scaffoldStrainL:
        for strain2 in scaffoldStrainL:
            blastFnL.append(os.path.join(blastDir,strain1+'_-VS-_'+strain2+'.out'))
    
    # calc scores for scaffold strains
    scoresO=Score.Score()
    scoresO.initializeDataAttributes(blastFnL,paramD,tuple(scaffoldStrainL))
    scoresO = scores.calcRawScores(paramD,scoresO)
    scoresO = scores.calcSynScores(scoresO,geneOrderD,paramD)
    scoresO = scores.calcCoreSynScores(scoresO,scaffoldStrainL,paramD,geneOrderD)    
    scores.writeScores(scoresO,scaffoldStrainL,paramD['scoresFN'])

    # families
    scaffoldAabrhL = scores.loadOrthos(paramD['aabrhFN'])
    with open(paramD['familyFormationSummaryFN'],'w') as familyFormationSummaryF:
        initialFamiliesO,originFamiliesO = families.createFamiliesO(scaffoldTree,scaffoldStrainL,scoresO,genesO,scaffoldAabrhL,paramD,familyFormationSummaryF)

    return initialFamiliesO,originFamiliesO

def blastRepGenesVsScaffold(paramD,originFamiliesO,scaffoldTree,scaffoldFamilyRepGenesFastaFN,scaffoldFamilyRepGenesNumPerFamily,genesO,allStrainsT,fastaDir):
    '''Get representative genes from each locus family, put in a db, and
blast vs. all the strains.
    '''
    
    seqD = genomes.loadSeq(paramD,'_prot.fa',originFamiliesO.getAllGenes())

    # get a set of representative genes from each locusfamily
    repGenesIterator = getGeneSubsetFromLocusFamilies(originFamiliesO, scaffoldTree, scaffoldFamilyRepGenesNumPerFamily, genesO)

    # write to db
    with open(scaffoldFamilyRepGenesFastaFN,"w") as f:
        for geneNum in repGenesIterator:
            writeFastaLine(f,int(geneNum),seqD)

    # Blast them vs. all strains
    allStrainsFileNamesL = []
    for strain in allStrainsT:
        allStrainsFileNamesL.append(os.path.join(fastaDir,strain+"_prot.fa"))
    
    blast.runBlast([scaffoldFamilyRepGenesFastaFN],allStrainsFileNamesL,paramD)

# old families.py functions
def getGeneSubsetFromLocusFamilies(originFamiliesO,tree,numRepresentativeGenesPerLocFam,genesO):
    '''Loop over all locus families and sample
numRepresentativeGenesPerLocFam from each.'''

    # get some stuff from tree
    leafL = tree.leaves()
    subtreeD= tree.createSubtreeD()
    numNodes = tree.nodeCount()
    
    for lfO in originFamiliesO.iterLocusFamilies():
        for geneNum in getGeneSubsetFromOneLocusFamily(lfO,originFamiliesO,numRepresentativeGenesPerLocFam,leafL,numNodes,genesO,subtreeD):
            yield geneNum

def getGeneSubsetFromOneLocusFamily(lfO,originFamiliesO,numRepresentativeGenesPerLocFam,leafL,numNodes,genesO,subtreeD):
    '''Sample numRepresentativeGenesPerLocFam genes from one
LocusFamily. This function simply divides the LocusFamily genes into
sets on the left and right branches, and attempts to take a similar
sized sample from each.'''

    lfGenesL = list(lfO.iterGenes())
    for i in range(len(lfGenesL)):
        lfGenesL[i] = str(lfGenesL[i])
    if len(lfGenesL) <= numRepresentativeGenesPerLocFam:
        return lfGenesL
    elif lfO.lfMrca in leafL:
        # it's a tip
        return random.sample(lfGenesL,numRepresentativeGenesPerLocFam)
    else:
        famO = originFamiliesO.getFamily(lfO.famNum)
        famTreeO = famO.geneTreeO
        if famTreeO != None:
            lfTreeO = famTreeO.prune(lfGenesL)
            return pruneGeneTree(lfTreeO, numRepresentativeGenesPerLocFam)
        else:
            return random.sample(lfGenesL,numRepresentativeGenesPerLocFam)
   
def mapAllStrainsGenes(paramD,numGenesAllStrains,originFamiliesO,allStrainsT,blastDir):
    '''For every gene in all the strains, map onto the scaffold genes
using the blast output. Creates an array, where the index is gene
number and the value is family number. unmappedVal is reserved to
indicate no blast similarity.
    '''
    # decide 32 or 16 bit for locFamNumAr
    numFam = originFamiliesO.getNumFamilies()
    if numFam > 2**16:
        print("Number of families is",numFam,"so using 32 bit array",file=sys.stderr)
        arSizeInBytes = 4
        unmappedVal = 2**32-1
        locFamNumAr = numpy.full(numGenesAllStrains, fill_value = unmappedVal, dtype=ctypes.c_uint32)
    else:
        arSizeInBytes = 2
        unmappedVal = 2**16-1
        locFamNumAr = numpy.full(numGenesAllStrains, fill_value = unmappedVal, dtype=ctypes.c_ushort)
        
    evalueAr = numpy.ones(numGenesAllStrains, dtype=ctypes.c_float)


    # first, go through all genes in scaffold strains and correctly
    # mark their family in the array. at same time create gene2FamNumD
    gene2LocFamNumD = {}
    for locusFamO in originFamiliesO.iterLocusFamilies():
        for gene in locusFamO.iterGenes():
            locFamNum = locusFamO.locusFamNum
            locFamNumAr[gene] = locFamNum
            evalueAr[gene] = 0 # so won't be changed.
            gene2LocFamNumD[int(gene)] = locFamNum
            
    # fill array with family values, indicating mapping for each gene
    for strain in allStrainsT:
        scaffoldFamilyRepGenesFastaStem = paramD['scaffoldFamilyRepGenesFastaFN'].split(".")[0]
        blastFN = os.path.join(blastDir,scaffoldFamilyRepGenesFastaStem+'_-VS-_'+strain+'.out')
        for queryGene,targetGene,evalue,alCov,pIdent,score in blast.parseBlastFile(blastFN,paramD['evalueThresh'],paramD['xlMapAlignCoverThresh'], paramD['percIdentThresh']):
            queryLocFam = gene2LocFamNumD[queryGene]
            if evalue < evalueAr[targetGene]:
                evalueAr[targetGene] = evalue
                locFamNumAr[targetGene] = queryLocFam

    # write array to file
    with open(paramD['allGenesMapToScaffoldFN'],'wb') as f:
        # initial 8 bytes indicate the array type (size in bits)
        f.write(arSizeInBytes.to_bytes(8,'little'))
        # next 8 bytes are the length of the array
        f.write(numGenesAllStrains.to_bytes(8,'little'))

        # write out array
        for i in range(numGenesAllStrains):
            f.write(int(locFamNumAr[i]).to_bytes(arSizeInBytes,'little'))
        
    return locFamNumAr,unmappedVal
    
def readAllGenesMapToScaffoldAr(allGenesMapToScaffoldFN):
    '''Read in binary array with mapping from all genes to scaffold.'''

    with open(allGenesMapToScaffoldFN,'rb') as f:

        arSizeInBytes = int.from_bytes(f.read(8),'little')
        numGenesAllStrains = int.from_bytes(f.read(8),'little')

        if arSizeInBytes == 4:
            unmappedVal = 2**32-1
            locFamNumAr = numpy.full(numGenesAllStrains, fill_value = unmappedVal, dtype=ctypes.c_uint32)
        else:
            unmappedVal = 2**16-1
            locFamNumAr = numpy.full(numGenesAllStrains, fill_value = unmappedVal, dtype=ctypes.c_ushort)

        # now fill values
        for i in range(numGenesAllStrains):
            locFamNumAr[i] = int.from_bytes(f.read(arSizeInBytes),'little')
           
    return locFamNumAr,unmappedVal

def writeUnmappedToFile(paramD,locFamNumAr,unmappedVal,fastaDir):
    '''Determine which genes are unmapped, and write them to unMappedGenesFN.'''
    # parameters
    numUnmappedGenesPerFasta = paramD['numUnmappedGenesPerFasta']
    unMappedGenesFilePathStem = paramD['unMappedGenesFilePathStem']
    
    # get unmapped
    unmappedGenesS = set()
    for gn in range(len(locFamNumAr)):
        if locFamNumAr[gn] == unmappedVal:
            # not mapped, save to db
            unmappedGenesS.add(gn)

    # load seqs
    seqD = genomes.loadSeq(paramD,'_prot.fa',unmappedGenesS)

    # get file name list for putting unmapped genes in
    numFiles = math.ceil(len(unmappedGenesS) / numUnmappedGenesPerFasta)
    fnL = [os.path.join(fastaDir,unMappedGenesFilePathStem+str(i)+'.fa') for i in range(numFiles)]

    # put unmapped genes in numFiles sublists
    unmappedGeneGroupsL = []
    subL = []
    numGenes = len(unmappedGenesS)
    for i in range(numGenes):
        gn = unmappedGenesS.pop()
        subL.append(gn)
        if len(subL) == numUnmappedGenesPerFasta:
            unmappedGeneGroupsL.append(subL)
            subL = []
    unmappedGeneGroupsL.append(subL) # add final one
            
    # write to file
    i=0
    for fn in fnL:
        with open(fn,'w') as f:
            for geneNum in unmappedGeneGroupsL[i]:
                writeFastaLine(f,geneNum,seqD)
        i+=1
        
    return fnL
                
def pickAdditionalStrains(allStrainsT,scaffoldStrainsL,blastDir,unMappedGenesFileL,evalueThresh,alignCoverThresh,genesO,numStrainsToAddToScaffold,percIdentThresh):
    '''Given blast output of unMapped genes against itself, pick strains
which will maximize the number of additional genes we can map.
    '''
    # create dict to hold counts of hits
    #hitCountD={}
    #for strain in allStrainsT:
    #    if strain not in scaffoldStrainsL:
    #        hitCountD[strain] = 0

    if numStrainsToAddToScaffold == 0:
        return []
            
    hitArray = [[0 for i in range(len(allStrainsT))] for x in range(len(allStrainsT))]

    # parse blast
    for unMappedGenesFN in unMappedGenesFileL:
        fileOnly = os.path.split(unMappedGenesFN)[-1]
        unMappedGenesFileStem = fileOnly.split(".")[0]
        blastFN = os.path.join(blastDir,unMappedGenesFileStem+'_-VS-_'+unMappedGenesFileStem+'.out')
        for queryGene,targetGene,evalue,alCov,pIdent,score in blast.parseBlastFile(blastFN,evalueThresh,alignCoverThresh,percIdentThresh):
            strain1=genesO.numToStrainName(queryGene)
            strainNum1 = allStrainsT.index(strain1)
            strain2=genesO.numToStrainName(targetGene)
            strainNum2 = allStrainsT.index(strain2)
            hitArray[strainNum1][strainNum2] += 1
            hitArray[strainNum2][strainNum1] += 1

    strainsToAddL = []
    while len(strainsToAddL) < numStrainsToAddToScaffold:
        maxStrain = 0
        maxSum = 0
        for i in range(len(hitArray)):
            total = sum(hitArray[i])
            if total > maxSum:
                maxSum = total
                maxStrain = i
        strainToAdd = allStrainsT[maxStrain]
        strainsToAddL.append(strainToAdd)
        for i in range(len(hitArray)):
            hitArray[i][maxStrain] = 0
        hitArray[maxStrain] = [0 for i in range(len(allStrainsT))]
    
    return strainsToAddL
        
###### Print analysis xl mode ######

def xlAnalysisSummary(paramD):
    '''Summarize the mapping from all genes to the scaffold.'''

    # file name
    xlAnalysisSummaryFN = os.path.join(paramD['analysisDir'],paramD['xlAnalysisSummaryStem']+'.txt')

    # read array
    locFamNumAr,unmappedVal = readAllGenesMapToScaffoldAr(paramD['allGenesMapToScaffoldFN'])

    # run over it
    mappedCt,unmappedCt,mappedProp,unmappedProp = arrayCounts(locFamNumAr,unmappedVal)
    
    # save to file
    rowsL=[]
    rowsL.append(["","Mapped genes","Unmapped genes"])
    rowsL.append(["Count",str(mappedCt),str(unmappedCt)])
    rowsL.append(["Proportion of whole",format(mappedProp,".3f"),format(unmappedProp,".3f")])

    with open(xlAnalysisSummaryFN,'w') as f:
        analysis.printTable(rowsL,indent=2,fileF=f)

def arrayCounts(locFamNumAr,unmappedVal):
    '''Count mapped and unmapped genes in the mapping array.'''

    numGenes = len(locFamNumAr)
    unmappedCt = 0
    for gn in range(numGenes):
        if locFamNumAr[gn] == unmappedVal:
            unmappedCt +=1

    mappedCt = numGenes-unmappedCt
    mappedProp = mappedCt/numGenes
    unmappedProp = unmappedCt/numGenes

    return mappedCt,unmappedCt,mappedProp,unmappedProp
