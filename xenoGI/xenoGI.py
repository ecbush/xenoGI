"""Provides the entry point to xenoGI's functionality."""
__version__ = "3.1.1"
import sys, glob, os, readline, rlcompleter
from . import parameters,genbank,blast,trees,genomes,Score,scores,families,islands,analysis,islandBed
from .Tree import *
from .Family import *

def main():
    
    #### check command line
    try:
        assert(len(sys.argv) == 3)
        paramFN=sys.argv[1]
        task = sys.argv[2]
        assert(task in ['parseGenbank', 'runBlast', 'calcScores','makeSpeciesTree', 'makeFamilies', 'makeIslands','refine', 'printAnalysis', 'interactiveAnalysis', 'createIslandBed', 'plotScoreHists', 'aminoAcidIdentity', 'runAll', 'version', 'debug'])
    
    except:
        print(
            """
   Exactly two arguments required.
      1. A path to a parameter file.
      2. The task to be run which must be one of: parseGenbank, runBlast, calcScores, makeSpeciesTree, makeFamilies, makeIslands, refine, printAnalysis, interactiveAnalysis, createIslandBed, plotScoreHists, aminoAcidIdentity, runAll or version.

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

    #### refineFamilies
    elif task == 'refine':
        refineWrapper(paramD)
        
    #### printAnalysis
    elif task == 'printAnalysis':
        printAnalysisWrapper(paramD,paramD['speciesTreeFN'],paramD['rootFocalClade'])

    #### interactiveAnalysis
    elif task == 'interactiveAnalysis':
        interactiveAnalysisWrapper(paramD)
        
    #### createIslandBed
    elif task == 'createIslandBed':
        createIslandBedWrapper(paramD)

    #### plotScoreHists
    elif task == 'plotScoreHists':
        scores.plotScoreHists(paramD)

    #### aminoAcidIdentity
    elif task == 'aminoAcidIdentity':
        '''Assumes runBlast has already been run.'''
        strainNamesT = readStrainInfoFN(paramD['strainInfoFN'])
        aaiD=analysis.calculateAAI(paramD,strainNamesT)
        analysis.printAminoAcidIdentity(aaiD,strainNamesT)
        
    #### runAll
    elif task == 'runAll':
        parseGenbankWrapper(paramD)
        runBlastWrapper(paramD)
        blastFnL=glob.glob(paramD['blastFilePath'])
        calcScoresWrapper(paramD,blastFnL)
        makeFamiliesWrapper(paramD)
        makeIslandsWrapper(paramD)
        refineWrapper(paramD)
        printAnalysisWrapper(paramD,paramD['speciesTreeFN'],paramD['rootFocalClade'])
        createIslandBedWrapper(paramD)

    #### version
    elif task == 'version':
        print("xenoGI",__version__)
        
    #### debug
    elif task == 'debug':
        debugWrapper(paramD)

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
    scoresO.initializeDataAttributes(blastFnL,paramD,strainNamesT,genesO)

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

def loadTreeRelatedData(speciesTreeFN):
    """Load some data related to trees."""
    speciesRtreeO = Rtree()
    speciesRtreeO.fromNewickFileLoadSpeciesTree(speciesTreeFN)
    subtreeD=speciesRtreeO.createSubtreeD()
    return speciesRtreeO,subtreeD

def makeFamiliesWrapper(paramD):
    """Wrapper to create gene families."""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    speciesRtreeO,subtreeD = loadTreeRelatedData(paramD['speciesTreeFN'])

    ## read scores
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    aabrhHardCoreL = scores.loadOrthos(paramD['aabrhFN'])

    ## make gene families
    with open(paramD['familyFormationSummaryFN'],'w') as familyFormationSummaryF:
        initialFamiliesO,originFamiliesO = families.createFamiliesO(speciesRtreeO,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,familyFormationSummaryF)
        
def makeIslandsWrapper(paramD):
    """Wrapper to create islands"""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    speciesRtreeO,subtreeD = loadTreeRelatedData(paramD['speciesTreeFN'])
    
    ## read gene families
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],speciesRtreeO,genesO,"origin")

    ## group gene families into islands
    with open(paramD['islandFormationSummaryFN'],'w') as islandFormationSummaryF:
        locIslByNodeD = islands.makeLocusIslands(geneOrderD,subtreeD,speciesRtreeO,paramD,originFamiliesO,paramD['rootFocalClade'],islandFormationSummaryF)

def refineWrapper(paramD):
    """Wrapper to refine families and islands by considering alternate reconcilitions."""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    speciesRtreeO,subtreeD = loadTreeRelatedData(paramD['speciesTreeFN'])
    
    ## read gene families
    initialFamiliesO = families.readFamilies(paramD['initFamilyFN'],speciesRtreeO,genesO,"initial")
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],speciesRtreeO,genesO,"origin")
    islandByNodeD=islands.readIslands(paramD['islandOutFN'],speciesRtreeO)

    ## refine families
    with open(paramD['familyFormationSummaryFN'],'a') as familyFormationSummaryF:
        initialFamiliesO,originFamiliesO = families.refineFamilies(paramD,islandByNodeD,initialFamiliesO,originFamiliesO,geneOrderD,genesO,familyFormationSummaryF,strainNamesT)

    ## redo islands with new originFamilies
    with open(paramD['islandFormationSummaryFN'],'w') as islandFormationSummaryF:
        locIslByNodeD = islands.makeLocusIslands(geneOrderD,subtreeD,speciesRtreeO,paramD,originFamiliesO,paramD['rootFocalClade'],islandFormationSummaryF)
        
def printAnalysisWrapper(paramD,speciesTreeFN,rootFocalClade):
    """Wrapper to run analysis. We take speciesTreeFN and rootFocalClade as
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
    speciesRtreeO,subtreeD = loadTreeRelatedData(speciesTreeFN)
    strainNamesT=speciesRtreeO.leaves()
    geneOrderD=genomes.createGeneOrderD(paramD['geneOrderFN'],strainNamesT)
    genesO = genomes.genes(paramD['geneInfoFN'])
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],speciesRtreeO,genesO,"origin")
    islandByNodeD=islands.readIslands(paramD['islandOutFN'],speciesRtreeO)
    
    ## analysis

    # Print islands summary file
    with open(islandsSummaryFN,'w') as islandsOutF:
        analysis.vPrintAllLocusIslands(islandByNodeD,speciesRtreeO,rootFocalClade,subtreeD,originFamiliesO,genesO,islandsOutF)

    # print islands tsv file
    with open(islandsTsvFN,'w') as islandsTsvF:
        analysis.printAllLocusIslandsTsv(islandByNodeD,speciesRtreeO,rootFocalClade,originFamiliesO,genesO,islandsTsvF)
    
    # Print species files with all the genes, grouped by contig
    gene2FamIslandD = analysis.createGene2FamIslandD(islandByNodeD,originFamiliesO)
    analysis.printSpeciesContigs(geneOrderD,genesFNstem,".tsv",genesO,gene2FamIslandD,originFamiliesO,rootFocalClade,strainNamesT)

def createIslandBedWrapper(paramD):
    """Wrapper to make output bed files."""

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    speciesRtreeO,subtreeD = loadTreeRelatedData(paramD['speciesTreeFN'])
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],speciesRtreeO,genesO,"origin")
    islandByNodeD = islands.readIslands(paramD['islandOutFN'],speciesRtreeO)
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    genesO.initializeGeneNumToNameD(paramD['geneInfoFN'],strainNamesT)
    
    # get islands organized by strain
    islandByStrainD = islandBed.createIslandByStrainD(strainNamesT,islandByNodeD,originFamiliesO,genesO)

    islandBed.createAllBeds(islandByStrainD,genesO,speciesRtreeO,strainNamesT,paramD,originFamiliesO)

def interactiveAnalysisWrapper(paramD):
    """Enter interactive mode."""

    ## Set up the modules a bit differently for interactive mode
    import code,sys
    from xenoGI.analysis import createGene2FamIslandD,printScoreMatrix,matchFamilyIsland,printLocusIslandNeighb,vPrintLocusIslandsAtNode,printOutsideFamilyScores,getOFamsWithEvent,countAllEvents
    from .xenoGI import families
    
    ## Wrapper analysis functions. For convenience these assume a
    ## bunch of global variables.
    
    def printFam(familiesO,familyNum,fileF=sys.stdout):
        '''This is a wrapper to provide an easy way to print relevant info on
    a family. For ease of use, we take only two arguments, assuming
    all the other required stuff is available at the top
    level. familyNum is the numerical identifier of a family, and
    familiesO is a family object.
        '''

        print("Family",familyNum,file=fileF)
        fam = familiesO.getFamily(familyNum)
        # print out the locus families
        for lfO in fam.getLocusFamilies():
            print("    LocusFamily",lfO.locusFamNum,file=fileF)
            print("      lfMrca "+lfO.lfMrca,file=fileF)
            for geneNum in lfO.iterGenes():
                print("      "+genesO.numToName(geneNum),file=fileF)
            
        if hasattr(fam,"sourceFam"):
            print(file=fileF)
            print("    Source family",fam.sourceFam,file=fileF)
            print(file=fileF)
            
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

        if hasattr(fam,"printReconByGeneTree"):
            # is ofam
            ofam = fam
        else:
            # is ifam, we must convert dtlorMprD to node format before printing
            speciesPreOrderT = speciesRtreeO.preorder()
            dtlorMprD = fam.getMprReconDFromMpr(speciesPreOrderT,paramD)
            ofam = originFamily(None,None,geneTreeO=fam.geneTreeO,dtlorMprD=dtlorMprD)
            ofam.printReconByGeneTree(genesO)

        print("Gene tree",file=fileF)
        print(ofam.geneTreeO.toNewickStr(),file=fileF)
        print()
        print("Gene tree annotated with reconciliation [branch events | node events]",file=fileF)
        print(ofam.getNewickGeneTreeWithReconLabels(genesO,includeBrLength=True),file=fileF)
        print(file=fileF)
        print("Reconciliation of gene tree onto species tree",file=fileF)       
        ofam.printReconByGeneTree(genesO,fileF)
            
        
    def findGene(searchStr,fileF=sys.stdout):
        '''Find information about a gene. Searches all the fields present in
    the geneInfo file, so the search string can be a locus tag, protein ID, a
    common name, or something present in the description. For each
    hit, prints the gene, LocusIsland, initialFamily, originFamily,
    LocusFamily and gene description. This is a wrapper that assumes
    various required objects are present at the top level.
        '''
        L=matchFamilyIsland(genesO,gene2FamIslandD,searchStr)
        for geneName,locusIslandNum,ifamNum,ofamNum,locusFamNum,descrip in L:
            print("<gene:"+str(geneName),"locIsl:"+str(locusIslandNum),"ifam:"+str(ifamNum),"ofam:"+str(ofamNum),"locFam:"+str(locusFamNum),descrip+">",file=fileF)

    def printLocusIsland(locusIslandNum,synWSize,fileF=sys.stdout):
        '''Print a LocusIsland and its genomic context in each species. We
        include synWSize/2 genes in either direction beyond the locus island.
        '''
        printLocusIslandNeighb(locusIslandNum,synWSize,subtreeD,islandByNodeD,originFamiliesO,geneOrderD,gene2FamIslandD,genesO,paramD['rootFocalClade'],fileF)


    def printLocusIslandsAtNode(node,fileF=sys.stdout):
        '''This is a wrapper to provide an easy way to print all the LocusIslands
    at a particular node in the speciesRtreeO. For ease of use, we take only a node
    number as argument, assuming all the other required stuff is available
    at the top level.
        '''
        vPrintLocusIslandsAtNode(islandByNodeD[node],paramD['rootFocalClade'],subtreeD,originFamiliesO,genesO,fileF)

    ## Load data

    print("Loading data, this may take some time...")
    
    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    speciesRtreeO,subtreeD = loadTreeRelatedData(paramD['speciesTreeFN'])
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    genesO.initializeGeneNumToNameD(paramD['geneInfoFN'],strainNamesT)
    initialFamiliesO = families.readFamilies(paramD['initFamilyFN'],speciesRtreeO,genesO,"initial")
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],speciesRtreeO,genesO,"origin")
    islandByNodeD=islands.readIslands(paramD['islandOutFN'],speciesRtreeO)
    gene2FamIslandD = createGene2FamIslandD(islandByNodeD,originFamiliesO)
    scoresO = scores.readScores(strainNamesT,paramD['scoresFN'])
    scoresO.createNodeConnectD() # make nodeConnectD attribute
    
    # set up interactive console
    variables = globals()
    variables.update(locals())
    readline.set_completer(rlcompleter.Completer(variables).complete)
    readline.parse_and_bind("tab: complete")
    print("Ready for interactive analysis.")
    code.InteractiveConsole(variables).interact()

def debugWrapper(paramD):
    '''Take us into interactive mode for debugging.'''

    ## Set up the modules a bit differently for interactive mode
    import code,sys,numpy
    from .xenoGI import analysis,parameters,trees,genomes,families,islands,analysis,Score,scores

    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    speciesRtreeO,subtreeD = loadTreeRelatedData(paramD['speciesTreeFN'])
    
    # set up interactive console
    variables = globals()
    variables.update(locals())
    readline.set_completer(rlcompleter.Completer(variables).complete)
    readline.parse_and_bind("tab: complete")
    code.InteractiveConsole(variables).interact()
