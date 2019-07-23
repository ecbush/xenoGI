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
        assert(task in ['parseGenbank', 'runBlast', 'calcScores','makeSpeciesTree', 'makeFamilies', 'makeIslands', 'printAnalysis', 'createIslandBed', 'plotScoreHists', 'interactiveAnalysis', 'runAll','makeGeneFamilyTrees', 'version', 'debug', 'simValidation'])
    
    except:
        print(
            """
   Exactly two arguments required.
      1. A path to a parameter file.
      2. The task to be run which must be one of: parseGenbank, runBlast, calcScores, makeSpeciesTree, makeFamilies, makeIslands, printAnalysis, createIslandBed, plotScoreHists, interactiveAnalysis, runAll, makeGeneFamilyTrees or version.

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
        printAnalysisWrapper(paramD)

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
        printAnalysisWrapper(paramD)
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
    scoresO.initializeDataAttributes(blastFnL)

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

def loadTreeRelatedData(paramD):
    """Load some data related to trees."""
    tree = trees.readTree(paramD['treeFN'])
    subtreeD=trees.createSubtreeD(tree)
    return tree,subtreeD

def makeFamiliesWrapper(paramD):
    """Wrapper to create gene families."""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD)

    ## read scores
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    aabrhHardCoreL = scores.loadOrthos(paramD['aabrhFN'])

    ## make gene families
    with open(paramD['familyFormationSummaryFN'],'w') as familyFormationSummaryF:
        familiesO = families.createFamiliesO(tree,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,subtreeD,familyFormationSummaryF)

def makeIslandsWrapper(paramD):
    """Wrapper to create islands"""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD)
    
    ## read gene families
    familiesO = families.readFamilies(paramD['familyFN'],tree,genesO)

    ## group gene families into islands
    with open(paramD['islandFormationSummaryFN'],'w') as islandFormationSummaryF:
        locIslByNodeD = islands.makeLocusIslands(geneOrderD,subtreeD,tree,paramD,familiesO,islandFormationSummaryF)

def printAnalysisWrapper(paramD):
    """Wrapper to run analysis."""

    ## setup output directory and file names
    # if directory for analysis doesn't exist yet, make it
    analDir = paramD['analysisFilePath'].split("*")[0]
    if glob.glob(analDir)==[]:
        os.mkdir(analDir)

    islandSummaryStem = paramD['islandsSummaryStem']
    analExtension = paramD['analysisFilePath'].split("*")[1]
    islandsSummaryFN = os.path.join(analDir,islandSummaryStem+analExtension)
    genesFNstem = os.path.join(analDir,paramD['genesFNstem'])
    
    ## load stuff
    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD)
    
    islandByNodeD=islands.readIslands(paramD['islandOutFN'],tree)
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    familiesO = families.readFamilies(paramD['familyFN'],tree,genesO)

    ## analysis

    # Print out all islands
    islandsOutF = open(islandsSummaryFN,'w')
    analysis.vPrintAllLocusIslands(islandByNodeD,tree,paramD['rootFocalClade'],subtreeD,familiesO,genesO,islandsOutF)
    islandsOutF.close()

    # Print species files with all the genes, grouped by contig
    gene2FamIslandD = analysis.createGene2FamIslandD(islandByNodeD,familiesO)
    analysis.printSpeciesContigs(geneOrderD,genesFNstem,analExtension,genesO,gene2FamIslandD,familiesO,strainNamesT)

def createIslandBedWrapper(paramD):
    """Wrapper to make output bed files."""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD)
    familiesO = families.readFamilies(paramD['familyFN'],tree,genesO)
    islandByNodeD = islands.readIslands(paramD['islandOutFN'],tree)
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    genesO.initializeGeneNumToNameD(paramD['geneInfoFN'],strainNamesT)
    
    # get islands organized by strain
    islandByStrainD = islandBed.createIslandByStrainD(strainNamesT,islandByNodeD,familiesO,genesO)

    islandBed.createAllBeds(islandByStrainD,genesO,tree,strainNamesT,paramD)

def interactiveAnalysisWrapper(paramD):
    """Enter interactive mode."""

    ## Set up the modules a bit differently for interactive mode
    import code,sys
    from xenoGI.analysis import createGene2FamIslandD,printScoreMatrix,matchFamilyIsland,printLocusIslandNeighb,vPrintLocusIslandsAtNode,coreNonCoreCtAtNode,printOutsideFamilyScores

    ## Wrapper analysis functions. For convenience these assume a
    ## bunch of global variables.
    
    def printFam(familyNum,fileF=sys.stdout):
        '''This is a wrapper to provide an easy way to print relevant info on
    a family. For ease of use, we take only one argument, assuming all the
    other required stuff is available at the top level. familyNum is the
    numerical identifier of a family.
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
        printLocusIslandNeighb(locusIslandNum,synWSize,subtreeD,islandByNodeD,familiesO,geneOrderD,gene2FamIslandD,genesO,fileF)


    def printLocusIslandsAtNode(node,fileF=sys.stdout):
        '''This is a wrapper to provide an easy way to print all the LocusIslands
    at a particular node in the tree. For ease of use, we take only a node
    number as argument, assuming all the other required stuff is available
    at the top level.
        '''
        vPrintLocusIslandsAtNode(islandByNodeD[node],subtreeD,familiesO,genesO,fileF)


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
            nonCore,core=coreNonCoreCtAtNode(tree,node,familyByNodeL,familiesO)
            rowL.append([node,str(core),str(nonCore),str(core+nonCore),str(format(nonCore/(core+nonCore),".3f"))])
        printTable(rowL,fileF=fileF)
        print(file=fileF)
        return


    ## Load data
    
    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD)    
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    genesO.initializeGeneNumToNameD(paramD['geneInfoFN'],strainNamesT)
    familiesO = families.readFamilies(paramD['familyFN'],tree,genesO)
    islandByNodeD=islands.readIslands(paramD['islandOutFN'],tree)
    gene2FamIslandD = createGene2FamIslandD(islandByNodeD,familiesO)
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    scoresO.createNodeConnectD() # make nodeConnectD attribute
    
    code.interact(local=locals())

def makeGeneFamilyTreesWrapper(paramD):
    '''Create a gene tree for each family.'''

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD)
    familiesO = families.readFamilies(paramD['familyFN'],tree,genesO)
    trees.makeGeneFamilyTrees(paramD,genesO,familiesO)
    
def debugWrapper(paramD):
    '''Take us into interactive mode for debugging.'''

    ## Set up the modules a bit differently for interactive mode
    import code,sys,numpy
    from .xenoGI import parameters,trees,genomes,families,islands,analysis,Score,scores

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    tree,subtreeD = loadTreeRelatedData(paramD)
    
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

