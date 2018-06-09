"""Provides the entry point to xenoGI's functionality."""
__version__ = "1.1.0"
import sys, glob, os
from . import parameters,genbank,blast,trees,genomes,Score,scores,Family,families,islands,analysis,islandBed

def main():
    
    #### check command line
    try:
        assert(len(sys.argv) == 3)
        paramFN=sys.argv[1]
        task = sys.argv[2]
        assert(task in ['parseGenbank', 'runBlast', 'calcScores', 'makeFamilies', 'makeIslands', 'printAnalysis', 'createIslandBed', 'runAll'])
    
    except:
        print(
            """
   Exactly two arguments required.
      1. A path to a parameter file.
            2. The task to be run which must be one of: parseGenbank, runBlast, calcScores, makeFamilies, makeIslands, printAnalysis, createIslandBed, or runAll.

   For example: 
      xenoGI params.py parseGenbank
"""
            ,file=sys.stderr)
        sys.exit(1)

        
    #### load parameters and some other data we'll use below
    paramD = parameters.loadParametersD(paramFN)
        
    #### parseGenbank
    if task == 'parseGenbank':
        parseGenbankWrapper(paramD)
        
    #### runBlast
    elif task == 'runBlast':
        blast.runBlast(paramD['fastaFilePath'],paramD['blastFilePath'],paramD['blastExecutDirPath'],paramD['blastCLine'],paramD['numThreads'])

    #### calcScores
    elif task == 'calcScores':
        calcScoresWrapper(paramD)

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
        
    #### runAll
    elif task == 'runAll':
        parseGenbankWrapper(paramD)
        blast.runBlast(paramD['fastaFilePath'],paramD['blastFilePath'],paramD['blastExecutDirPath'],paramD['blastCLine'],paramD['numThreads'])
        calcScoresWrapper(paramD)
        makeFamiliesWrapper(paramD)
        makeIslandsWrapper(paramD)
        printAnalysisWrapper(paramD)
        createIslandBedWrapper(paramD)

        
######## Task related functions

def parseGenbankWrapper(paramD):
    """Wrapper running stuff to parse genbank files."""
    genbankFileList=glob.glob(paramD['genbankFilePath'])

    # if directory for fastas doesn't exist yet, make it
    fastaDir = paramD['fastaFilePath'].split('*')[0]
    if glob.glob(fastaDir)==[]:
        os.mkdir(fastaDir)

    fileNameMapD = parameters.loadFileNameMapD(paramD['fileNameMapFN'],genbankFileList)

    # parse
    genbank.parseGenbank(paramD['geneOrderFN'],paramD['redundProtsFN'],paramD['geneInfoFN'],fastaDir,genbankFileList,fileNameMapD)

def loadMiscDataStructures(paramD):
    """Creates a few data structures that are used in multiple tasks."""

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    
    # an object for gene name conversions
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    return tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT



def calcScoresWrapper(paramD):
    """Wrapper running stuff to calculate scores."""

    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)

    # object for storing scores
    scoresO=Score.Score()
    scoresO.initializeDataAttributes(paramD['blastFilePath'],geneNames)

    ## similarity scores
    scoresO = scores.calcRawScores(paramD['fastaFilePath'],paramD['numThreads'],geneNames,paramD['gapOpen'],paramD['gapExtend'],paramD['matrix'],scoresO)

    ## normalized scores
    scoresO,aabrhL,aabrhRawScoreSummmaryD=scores.calcNormScores(tree,strainNum2StrD,paramD['blastFilePath'],paramD['evalueThresh'],scoresO,geneNames,paramD['aabrhFN'])

    ## synteny scores
    scoresO = scores.calcSynScores(scoresO,aabrhRawScoreSummmaryD,geneNames,geneOrderT,paramD['synWSize'],paramD['numSynToTake'],paramD['numThreads'])

    ## core synteny scores
    scoresO = scores.calcCoreSynScores(scoresO,aabrhL,geneNames,geneOrderT,paramD['coreSynWsize'])

    # write scores to file
    scores.writeScores(scoresO,geneNames,paramD['scoresFN'])
    
def makeFamiliesWrapper(paramD):
    """Wrapper to create gene families."""

    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)
    
    ## read scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)

    ## make gene families
    familyFormationSummaryF = open(paramD['familyFormationSummaryFN'],'w')
    familyL = families.families(tree,subtreeL,geneNames,scoresO,paramD['minNormThresh'],paramD['minCoreSynThresh'],paramD['minSynThresh'],paramD['synAdjustThresh'],paramD['synAdjustExtent'],paramD['familyFN'],strainNum2StrD,familyFormationSummaryF)
    familyFormationSummaryF.close()

def makeIslandsWrapper(paramD):
    """Wrapper to create islands"""

    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)
    
    ## read scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)

    ## read gene families
    familyL = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    ## group gene families into islands
    islandFormationSummaryF = open(paramD['islandFormationSummaryFN'],'w')
    islands.makeIslands(geneOrderT,geneNames,subtreeL,tree,paramD['proxThreshL'],familyL,paramD['numThreads'],strainStr2NumD,strainNum2StrD,paramD['rootFocalClade'],paramD['islandOutFN'],islandFormationSummaryF)
    islandFormationSummaryF.close()

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
    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)
    islandByNodeL=islands.readIslands(paramD['islandOutFN'],tree,strainStr2NumD)

    nodesL=trees.nodeList(tree)
    
    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])

    familyL = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    gene2FamD=analysis.createGene2FamD(familyL)
    fam2IslandD=analysis.createFam2IslandD(islandByNodeL)

    # scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)
    scoresO.createNodeConnectL(geneNames) # make nodeConnectL attribute

    ## analysis
    
    # calc family error scores
    families.calcErrorScores(familyL,scoresO,paramD['minNormThresh'],paramD['minCoreSynThresh'],paramD['minSynThresh'],paramD['famErrorScoreIncrementD'])

    # Print out all islands
    islandsOutF = open(islandsSummaryFN,'w')
    analysis.vPrintAllIslands(islandByNodeL,tree,paramD['rootFocalClade'],subtreeL,familyL,strainStr2NumD,strainNum2StrD,geneNames,geneInfoD,islandsOutF)
    islandsOutF.close()

    # Print species files with all the genes, grouped by contig
    analysis.printSpeciesContigs(geneOrderT,genesFNstem,analExtension,geneNames,gene2FamD,fam2IslandD,geneInfoD,familyL,strainNum2StrD)

def createIslandBedWrapper(paramD):
    """Wrapper to make output bed files."""

    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)
    
    leafNodesL = trees.leafList(tree)
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)
    familyL = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)
    islandByNodeL=islands.readIslands(paramD['islandOutFN'],tree,strainStr2NumD)
    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])    


    # get islands organized by strain
    islandByStrainD = islandBed.createIslandByStrainD(leafNodesL,strainNum2StrD,islandByNodeL,familyL,geneNames,geneInfoD)

    islandBed.createAllBeds(islandByStrainD,geneInfoD,tree,strainNum2StrD,paramD['bedFilePath'],paramD['scoreNodeMapD'],paramD['potentialRgbL'],paramD['bedNumTries'])

