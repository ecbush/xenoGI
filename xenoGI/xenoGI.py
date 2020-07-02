"""Provides the entry point to xenoGI's functionality."""
__version__ = "2.2.0"
import sys, glob, os
from . import parameters,genbank,blast,trees,genomes,Score,scores,Family,families,islands,analysis,islandBed

def main():
    
    #### check command line
    try:
        assert(len(sys.argv) == 3)
        paramFN=sys.argv[1]
        task = sys.argv[2]
        assert(task in ['parseGenbank', 'runBlast', 'calcScores','makeSpeciesTree', 'makeFamilies', 'makeFamiliesDTLOR', 'makeIslands', 'printAnalysis', 'createIslandBed', 'plotScoreHists', 'interactiveAnalysis', 'runAll','makeGeneFamilyTrees', 'version', 'debug', 'simValidation'])
    
    except:
        print(
            """
   Exactly two arguments required.
      1. A path to a parameter file.
      2. The task to be run which must be one of: parseGenbank, runBlast, calcScores, makeSpeciesTree, makeFamilies, makeFamiliesDTLOR, makeIslands, printAnalysis, createIslandBed, plotScoreHists, interactiveAnalysis, runAll, makeGeneFamilyTrees or version.

   For example: 
      xenoGI params.py parseGenbank
"""
            ,file=sys.stderr)
        sys.exit(1)

        
    #### load parameters and some other data we'll use below
    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)
        
    #### parseGenbank
    if task == 'parseGenbank':
        parseGenbankWrapper(paramD)
        
    #### runBlast
    elif task == 'runBlast':
        runBlastWrapper(paramD)
        
    #### calcScores
    elif task == 'calcScores':
        blastFnL=glob.glob(paramD['blastFilePath'])
        calcScoresWrapper(paramD,blastFnL)

    #### makeSpeciesTreeWrapper
    elif task == 'makeSpeciesTree':
        makeSpeciesTreeWrapper(paramD)
        
    #### makeFamilies
    elif task == 'makeFamilies':
        makeFamiliesWrapper(paramD)
        
    #### makeIslands
    elif task == 'makeIslands':
        makeIslandsWrapper(paramD)

    #### printAnalysis
    elif task == 'printAnalysis':
        printAnalysisWrapper(paramD,paramD['treeFN'],paramD['rootFocalClade'])

    #### createIslandBed
    elif task == 'createIslandBed':
        createIslandBedWrapper(paramD)

    #### interactiveAnalysis
    elif task == 'interactiveAnalysis':
        interactiveAnalysisWrapper(paramD)
        
    #### plotScoreHists
    elif task == 'plotScoreHists':
        scores.plotScoreHists(paramD)
        
    #### runAll
    elif task == 'runAll':
        parseGenbankWrapper(paramD)
        runBlastWrapper(paramD)
        blastFnL=glob.glob(paramD['blastFilePath'])
        calcScoresWrapper(paramD,blastFnL)
        makeFamiliesWrapper(paramD)
        makeIslandsWrapper(paramD)
        printAnalysisWrapper(paramD,paramD['treeFN'],paramD['rootFocalClade'])
        createIslandBedWrapper(paramD)

    #### makeGeneFamilyTrees
    elif task == 'makeGeneFamilyTrees':
        makeGeneFamilyTreesWrapper(paramD)
        
    #### version
    elif task == 'version':
        print("xenoGI",__version__)
        
    #### debug
    elif task == 'debug':
        debugWrapper(paramD)

    #### validation via simulation
    elif task == 'simValidation':
        simValidationWrapper(paramD)
        
######## Task related functions

def parseGenbankWrapper(paramD):
    """Wrapper running stuff to parse genbank files."""
    genbankFileList=glob.glob(paramD['genbankFilePath'])

    # if directory for fastas doesn't exist yet, make it
    fastaDir = paramD['fastaFilePath'].split('*')[0]
    if glob.glob(fastaDir)==[]:
        os.mkdir(fastaDir)

    fileNameMapD,strainNamesT = parameters.loadFileNameMapD(paramD['fileNameMapFN'],genbankFileList)
    writeStrainInfoFN(strainNamesT,paramD)
    
    # parse
    genbank.parseGenbank(paramD,fastaDir,genbankFileList,fileNameMapD)

def writeStrainInfoFN(strainNamesT,paramD):
    """Write strain numbers and names to strainInfoFN. Numbers are given
by the index in strainNamesT."""
    with open(paramD['strainInfoFN'],'w') as f:
        for strainNum in range(len(strainNamesT)):
            f.write(str(strainNum)+"\t"+strainNamesT[strainNum]+"\n")

def readStrainInfoFN(strainInfoFN):
    """Read strainNamesT from strainInfoFN."""
    strainNamesL=[]
    with open(strainInfoFN,'r') as f:
        while True:
            s=f.readline()
            if s=='':
                break
            strainNum,strainName = s.rstrip().split("\t")
            strainNamesL.append(strainName)
    return tuple(strainNamesL)
    
def runBlastWrapper(paramD):
    """Wrapper to blast all genome files against each other."""
    strainNamesT = readStrainInfoFN(paramD['strainInfoFN'])
    dbFileL=blast.getDbFileL(paramD['fastaFilePath'],strainNamesT)
    blast.runBlast(dbFileL,dbFileL,paramD)
    
def loadGenomeRelatedData(paramD):
    """Load some data related to genomes and strains."""
    strainNamesT = readStrainInfoFN(paramD['strainInfoFN'])
    genesO = genomes.genes(paramD['geneInfoFN'])
    geneOrderD=genomes.createGeneOrderD(paramD['geneOrderFN'],None)
    return strainNamesT,genesO,geneOrderD

def calcScoresWrapper(paramD,blastFnL):
    """Wrapper running stuff to calculate scores."""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    
    # object for storing scores
    scoresO=Score.Score()
    scoresO.initializeDataAttributes(blastFnL,paramD)

    ## similarity scores
    scoresO = scores.calcRawScores(paramD,scoresO)

    ## synteny scores
    scoresO = scores.calcSynScores(scoresO,geneOrderD,paramD)

    ## core synteny scores
    scoresO = scores.calcCoreSynScores(scoresO,strainNamesT,paramD,geneOrderD)

    # write scores to file
    scores.writeScores(scoresO,strainNamesT,paramD['scoresFN'])
    
def makeSpeciesTreeWrapper(paramD):
    '''call makeTree to create the species tree '''

    # need genesO should initializeGeneInfoD
    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    aabrhHardCoreL = scores.loadOrthos(paramD['aabrhFN'])
    trees.makeSpeciesTree(paramD,aabrhHardCoreL,genesO)

def loadTreeRelatedData(treeFN):
    """Load some data related to trees."""
    tree = trees.readTree(treeFN)
    subtreeD=trees.createSubtreeD(tree)
    return tree,subtreeD

def makeFamiliesWrapper(paramD):
    """Wrapper to create gene families."""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD['treeFN'])

    ## read scores
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    aabrhHardCoreL = scores.loadOrthos(paramD['aabrhFN'])

    ## make gene families
    with open(paramD['familyFormationSummaryFN'],'w') as familyFormationSummaryF:
        originFamiliesO = families.createFamiliesO(tree,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,familyFormationSummaryF)
        
def makeIslandsWrapper(paramD):
    """Wrapper to create islands"""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD['treeFN'])
    
    ## read gene families
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],tree,genesO)

    ## group gene families into islands
    with open(paramD['islandFormationSummaryFN'],'w') as islandFormationSummaryF:
        locIslByNodeD = islands.makeLocusIslands(geneOrderD,subtreeD,tree,paramD,originFamiliesO,paramD['rootFocalClade'],islandFormationSummaryF)

def printAnalysisWrapper(paramD,treeFN,rootFocalClade):
    """Wrapper to run analysis. We take treeFN and rootFocalClade as
arguments so we can pass in different things in different contexts
(e.g. in xlMode).
    """
    ## setup output directory and file names
    # if directory for analysis doesn't exist yet, make it
    analDir = paramD['analysisDir']
    if glob.glob(analDir)==[]:
        os.mkdir(analDir)

    islandFNStem = paramD['islandsFNStem']
    islandsSummaryFN = os.path.join(analDir,islandFNStem+"Summary.txt")
    islandsTsvFN = os.path.join(analDir,islandFNStem+".tsv")
    genesFNstem = os.path.join(analDir,paramD['genesFNstem'])
    
    ## load stuff
    tree,subtreeD = loadTreeRelatedData(treeFN)
    strainNamesT=tuple(trees.leafList(tree))
    geneOrderD=genomes.createGeneOrderD(paramD['geneOrderFN'],strainNamesT)
    genesO = genomes.genes(paramD['geneInfoFN'])
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],tree,genesO)
    islandByNodeD=islands.readIslands(paramD['islandOutFN'],tree)
    
    ## analysis

    # Print islands summary file
    with open(islandsSummaryFN,'w') as islandsOutF:
        analysis.vPrintAllLocusIslands(islandByNodeD,tree,rootFocalClade,subtreeD,originFamiliesO,genesO,islandsOutF)

    # print islands tsv file
    with open(islandsTsvFN,'w') as islandsTsvF:
        analysis.printAllLocusIslandsTsv(islandByNodeD,tree,rootFocalClade,originFamiliesO,genesO,islandsTsvF)
    
    # Print species files with all the genes, grouped by contig
    gene2FamIslandD = analysis.createGene2FamIslandD(islandByNodeD,originFamiliesO)
    analysis.printSpeciesContigs(geneOrderD,genesFNstem,".tsv",genesO,gene2FamIslandD,originFamiliesO,rootFocalClade,strainNamesT)

def createIslandBedWrapper(paramD):
    """Wrapper to make output bed files."""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD['treeFN'])
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],tree,genesO)
    islandByNodeD = islands.readIslands(paramD['islandOutFN'],tree)
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    genesO.initializeGeneNumToNameD(paramD['geneInfoFN'],strainNamesT)
    
    # get islands organized by strain
    islandByStrainD = islandBed.createIslandByStrainD(strainNamesT,islandByNodeD,originFamiliesO,genesO)

    islandBed.createAllBeds(islandByStrainD,genesO,tree,strainNamesT,paramD)

def interactiveAnalysisWrapper(paramD):
    """Enter interactive mode."""

    ## Set up the modules a bit differently for interactive mode
    import code,sys
    from xenoGI.analysis import createGene2FamIslandD,printScoreMatrix,matchFamilyIsland,printLocusIslandNeighb,vPrintLocusIslandsAtNode,coreNonCoreCtAtNode,printOutsideFamilyScores
    from .xenoGI import families # remove later
    
    ## Wrapper analysis functions. For convenience these assume a
    ## bunch of global variables.
    
    def printFam(familyNum,familiesO,fileF=sys.stdout):
        '''This is a wrapper to provide an easy way to print relevant info on
    a family. For ease of use, we take only two arguments, assuming
    all the other required stuff is available at the top
    level. familyNum is the numerical identifier of a family, and
    familiesO is a family object.
        '''
                    
        print("Family",familyNum,file=fileF)
        # print out the locus families
        for lfO in familiesO.getFamily(familyNum).getLocusFamilies():
            print("    LocusFamily",lfO.getStr(genesO," "),file=fileF)
        
        # print("Family error score (count of possibly misassigned genes):",familiesO[familyNum].possibleErrorCt,file=fileF)

        print(file=fileF)
        print("Matrix of raw similarity scores [0,1] between genes in the family",file=fileF)
        printScoreMatrix(familyNum,familiesO,genesO,scoresO,'rawSc',fileF)
        print(file=fileF)
        print(file=fileF)

        print("Matrix of core synteny scores [0,1] between genes in the family",file=fileF)
        printScoreMatrix(familyNum,familiesO,genesO,scoresO,'coreSynSc',fileF)
        print(file=fileF)
        print(file=fileF)

        print("Matrix of synteny scores [0,1] between genes in the family",file=fileF)
        printScoreMatrix(familyNum,familiesO,genesO,scoresO,'synSc',fileF)
        print(file=fileF)
        print(file=fileF)

        printOutsideFamilyScores(familyNum,familiesO,genesO,scoresO,fileF)
        print(file=fileF)
        print(file=fileF)

    def findLocusIsland(searchStr,fileF=sys.stdout):
        '''Print the gene, LocusIsland family and LocusFamily associated with
    searchStr. This is a wrapper that assumes various required objects
    are present at the top level.
        '''
        L=matchFamilyIsland(genesO,gene2FamIslandD,searchStr)
        for geneName,locusIslandNum, famNum, locusFamNum in L:
            print("<gene:"+str(geneName)+">","<locIsl:"+str(locusIslandNum)+">","<fam:"+str(famNum)+">","<locFam:"+str(locusFamNum)+">",file=fileF)

    def printLocusIsland(locusIslandNum,synWSize,fileF=sys.stdout):
        '''Print a LocusIsland and its genomic context in each species. We
        include synWSize/2 genes in either direction beyond the locus island.
        '''
        printLocusIslandNeighb(locusIslandNum,synWSize,subtreeD,islandByNodeD,originFamiliesO,geneOrderD,gene2FamIslandD,genesO,fileF)


    def printLocusIslandsAtNode(node,fileF=sys.stdout):
        '''This is a wrapper to provide an easy way to print all the LocusIslands
    at a particular node in the tree. For ease of use, we take only a node
    number as argument, assuming all the other required stuff is available
    at the top level.
        '''
        vPrintLocusIslandsAtNode(islandByNodeD[node],subtreeD,originFamiliesO,genesO,fileF)


    def printCoreNonCoreByNode(fileF=sys.stdout):
        '''For each node in the focal clade, print the number of core and
    non-core families. That is, for that node, we get all families present
    in descendent species. We then look at their mrcas. Those with an mrca
    below the node in question are non-core, others are core.'''

        focTree = trees.subtree(tree,paramD['rootFocalClade'])
        focInodesL=trees.iNodeList(focTree)
        familyByNodeL=createFamilyByNodeL(geneOrderD,gene2FamD)

        rowL=[]
        rowL.append(['Node','Core','Non-Core','Total','% Non-Core'])
        rowL.append(['----','----','--------','-----','----------'])
        for node in focInodesL:
            nonCore,core=coreNonCoreCtAtNode(tree,node,familyByNodeL,originFamiliesO)
            rowL.append([node,str(core),str(nonCore),str(core+nonCore),str(format(nonCore/(core+nonCore),".3f"))])
        printTable(rowL,fileF=fileF)
        print(file=fileF)
        return


    ## Load data
    
    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD['treeFN'])    
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    genesO.initializeGeneNumToNameD(paramD['geneInfoFN'],strainNamesT)
    initialFamiliesO = families.readFamilies(paramD['initFamilyFN'],tree,genesO)
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],tree,genesO)
    islandByNodeD=islands.readIslands(paramD['islandOutFN'],tree)
    gene2FamIslandD = createGene2FamIslandD(islandByNodeD,originFamiliesO)
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    scoresO.createNodeConnectD() # make nodeConnectD attribute
    
    code.interact(local=locals())

def makeGeneFamilyTreesWrapper(paramD):
    '''Create a gene tree for each family.'''

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD['treeFN'])
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],tree,genesO)
    trees.makeGeneFamilyTrees(paramD,genesO,originFamiliesO)
    
def debugWrapper(paramD):
    '''Take us into interactive mode for debugging.'''

    ## Set up the modules a bit differently for interactive mode
    import code,sys,numpy
    from .xenoGI import parameters,trees,genomes,families,islands,analysis,Score,scores
    from .xenoGI import Family
    
    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    speciesTree,subtreeD = loadTreeRelatedData(paramD['treeFN'])

    #initialFamiliesO = families.readFamilies(paramD['initFamilyFN'],speciesTree,genesO)
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],speciesTree,genesO)

    #ifam = initialFamiliesO.getFamily(1)
    
    D=float(paramD["duplicationCost"])
    T=float(paramD["transferCost"])
    L=float(paramD["lossCost"])
    O=float(paramD["originCost"])
    R=float(paramD["rearrangeCost"])
    rootFocalClade = paramD["rootFocalClade"]

    
    fam=originFamiliesO.getFamily(7506)
    fam.createGeneHistoryD(speciesTree,rootFocalClade)
    fam.geneHistoryD

    #iFamGeneTreeFileStem = 'initFam'
    #singleGeneInitFamNumL,multifurcatingL,bifurcatingL = families.loadGeneTrees(paramD,initialFamiliesO,iFamGeneTreeFileStem)
    #originFamiliesO = families.createOriginFamiliesO(speciesTree,singleGeneInitFamNumL,multifurcatingL,bifurcatingL,initialFamiliesO,genesO)
    
    """    
    for lfO in originFamiliesO.iterLocusFamilies():
        orig = lfO.origin(originFamiliesO)
        if orig == 'R':
            print(orig)
            lfO.printReconByGeneTree(originFamiliesO)
            print('---')


    with open("reconTemp.tsv","r") as f:
        s=f.readline() # skip first line
        while True:
            s=f.readline()
            if s=="":
                break
            initFamNum,optRootedGeneTree,optMPR,minCost,argT = s.rstrip().split("\t")
            initFamNum = eval(initFamNum)
            optMPR=eval(optMPR)
            optRootedGeneTree = eval(optRootedGeneTree)
            minCost=eval(minCost)
            argT=eval(argT)
            
            reconD = families.convertReconBranchToNode(optMPR,optRootedGeneTree)

            if initFamNum == 1:
                break


            ifam = initialFamiliesO.getFamily(initFamNum)
            #print(ifam.reconD)
            ifam.addGeneTree(optRootedGeneTree)
            ifam.addReconciliation(reconD)
            #print(ifam.reconD)
            #print("---")

    families.writeFamilies(initialFamiliesO,"initFamTemp.out",genesO,strainNamesT,paramD)
            

    print(families.costSum(reconD,D,T,L,O,R))

    families.parseOneBranchEntry(optMPR,('pTop', 'hTop', 4602),'*')

    for geneTree,splitReconD in families.getGeneTreeReconPairs(originTreeL,reconD):
        print(geneTree)
        print(splitReconD)
        print('---')
    """


    code.interact(local=locals())

    
def simValidationWrapper(paramD):
    '''Take us into interactive mode for running validation with
simulations. This is not really intended for the end user...'''

    ## Set up the modules a bit differently for interactive mode
    import sys
    from .xenoGI import parameters,trees,genomes,families,islands,analysis
    sys.path.insert(0,'/data/bushlab/htrans/dupDelProj/simCode')
    import validationSim

    moduleT=(parameters,trees,genomes,families,islands,analysis)
    
    xgiIslandGeneSetByNodeL,simIslandGeneSetByNodeL=validationSim.runValidation(moduleT,'simParams.py','params.py')
