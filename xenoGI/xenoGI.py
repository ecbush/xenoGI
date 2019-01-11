"""Provides the entry point to xenoGI's functionality."""
__version__ = "1.1.1"
import sys, glob, os
from . import parameters,genbank,blast,trees,genomes,Score,scores,Family,families,islands,analysis,islandBed

def main():
    
    #### check command line
    try:
        assert(len(sys.argv) == 3)
        paramFN=sys.argv[1]
        task = sys.argv[2]
        assert(task in ['parseGenbank', 'runBlast', 'calcScores', 'makeFamilies', 'makeIslands', 'printAnalysis', 'createIslandBed', 'plotScoreHists', 'interactiveAnalysis', 'runAll', 'debug'])
    
    except:
        print(
            """
   Exactly two arguments required.
      1. A path to a parameter file.
      2. The task to be run which must be one of: parseGenbank, runBlast, calcScores, makeFamilies, makeIslands, printAnalysis, createIslandBed, plotScoreHists, interactiveAnalysis, or runAll.

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
        blast.runBlast(paramD['fastaFilePath'],paramD['blastFilePath'],paramD['blastExecutDirPath'],paramD['blastCLine'],paramD['numThreads'],paramD['treeFN'])

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

    #### interactiveAnalysis
    elif task == 'interactiveAnalysis':
        interactiveAnalysisWrapper(paramD)
        
    #### plotScoreHists
    elif task == 'plotScoreHists':
        plotScoreHistsWrapper(paramD)

    #### runAll
    elif task == 'runAll':
        parseGenbankWrapper(paramD)
        blast.runBlast(paramD['fastaFilePath'],paramD['blastFilePath'],paramD['blastExecutDirPath'],paramD['blastCLine'],paramD['numThreads'],paramD['treeFN'])
        calcScoresWrapper(paramD)
        makeFamiliesWrapper(paramD)
        makeIslandsWrapper(paramD)
        printAnalysisWrapper(paramD)
        createIslandBedWrapper(paramD)

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
    scoresO.initializeDataAttributes(paramD['blastFilePath'],geneNames,strainStr2NumD)

    ## similarity scores
    scoresO = scores.calcRawScores(paramD['fastaFilePath'],paramD['numThreads'],geneNames,paramD['gapOpen'],paramD['gapExtend'],paramD['matrix'],scoresO)

    ## synteny scores
    scoresO = scores.calcSynScores(scoresO,geneNames,geneOrderT,paramD['synWSize'],tree,paramD['numSynToTake'],paramD['numThreads'])

    ## core synteny scores
    strainNamesL=sorted([strainNum2StrD[leaf] for leaf in trees.leafList(tree)])
    scoresO = scores.calcCoreSynScores(scoresO,strainNamesL,paramD['blastFilePath'],paramD['evalueThresh'],paramD['aabrhFN'],geneNames,geneOrderT,paramD['coreSynWsize'])
    
    # write scores to file
    scores.writeScores(scoresO,geneNames,paramD['scoresFN'])
    
def makeFamiliesWrapper(paramD):
    """Wrapper to create gene families."""

    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)
    ## read scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)
    aabrhL = scores.loadOrthos(paramD['aabrhFN'])

    ## make gene families
    familyFormationSummaryF = open(paramD['familyFormationSummaryFN'],'w')
    familiesO = families.createFamiliesO(tree,strainNum2StrD,scoresO,geneNames,aabrhL,paramD['singleStrainFamilyThresholdAdjust'],subtreeL,paramD['minCoreSynThresh'],paramD['synAdjustExtent'],familyFormationSummaryF,paramD['familyFN'])

    familyFormationSummaryF.close()

def makeIslandsWrapper(paramD):
    """Wrapper to create islands"""

    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)
    
    ## read scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)

    ## read gene families
    familiesO = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    ## group gene families into islands
    islandFormationSummaryF = open(paramD['islandFormationSummaryFN'],'w')
    locusIslandByNodeLMerged = islands.makeLocusIslands(geneOrderT,geneNames,subtreeL,tree,paramD['proxThreshL'],familiesO,paramD['numThreads'],strainStr2NumD,strainNum2StrD,paramD['rootFocalClade'],paramD['islandOutFN'],islandFormationSummaryF)

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

    scoresO = scores.readScores(paramD['scoresFN'],geneNames)

    ## analysis

    # Print out all islands
    islandsOutF = open(islandsSummaryFN,'w')
    analysis.vPrintAllIslands(islandByNodeL,tree,paramD['rootFocalClade'],subtreeL,familyL,strainStr2NumD,strainNum2StrD,geneNames,geneInfoD,islandsOutF)
    islandsOutF.close()

    # score calculations
    scoresO.createNodeConnectL(geneNames) # make nodeConnectL attribute
    families.calcErrorScores(familyL,scoresO,paramD['minNormThresh'],paramD['minCoreSynThresh'],paramD['minSynThresh'],paramD['famErrorScoreIncrementD'])

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

def plotScoreHistsWrapper(paramD):
    """Wrapper to make pdf of histograms of scores."""

    import matplotlib.pyplot as pyplot
    from matplotlib.backends.backend_pdf import PdfPages
    
    numBins = 80 # num bins in histograms
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)
    
    def scoreHists(scoresFN,outFN,strainNum2StrD,numBins,geneNames,scoreType):
        '''Read through a scores file, and separate into all pairwise comparisons. Then plot hist of each.'''

        # currently, this seems to require a display for interactive
        # plots. would be nice to make it run without that...

        scoresO = scores.readScores(scoresFN,geneNames)

        pyplot.ioff() # turn off interactive mode
        with PdfPages(outFN) as pdf:
            for strainPair in scoresO.getStrainPairs():
                fig = pyplot.figure()
                scoresL = list(scoresO.iterateScoreByStrainPair(strainPair,scoreType))
                pyplot.hist(scoresL,bins=numBins, density = True)
                pyplot.title(strainNum2StrD[strainPair[0]]+'-'+strainNum2StrD[strainPair[1]])
                pdf.savefig()
                pyplot.close()

    # plot histograms
    for scoreType,outFN in [('rawSc','rawSc.pdf'),('synSc','synSc.pdf'),('coreSynSc','coreSynSc.pdf'),]:
    
        scoreHists(paramD['scoresFN'],outFN,strainNum2StrD,numBins,geneNames,scoreType)

    
def interactiveAnalysisWrapper(paramD):
    """Enter interactive mode."""

    ## Set up the modules a bit differently for interactive mode
    import code,sys
    from xenoGI.analysis import createGene2FamD,createFam2IslandD,printScoreMatrix,matchFamilyIsland,printIslandNeighb,vPrintIslands,coreNonCoreCtAtNode,printOutsideFamilyScores,scoreHists,readScorePairs

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
            print("    LocusFamily",lfO.getStr(strainNum2StrD,geneNames," "),file=fileF)
        
        # print("Family error score (count of possibly misassigned genes):",familyL[familyNum].possibleErrorCt,file=fileF)

        print(file=fileF)
        print("Matrix of raw similarity scores [0,1] between genes in the family",file=fileF)
        printScoreMatrix(familyNum,subtreeL,familiesO,geneNames,scoresO,'rawSc',fileF)
        print(file=fileF)
        print(file=fileF)

        print(file=fileF)
        print("Matrix of normalized similarity scores between genes in the family",file=fileF)
        printScoreMatrix(familyNum,subtreeL,familiesO,geneNames,scoresO,'normSc',fileF)
        print(file=fileF)
        print(file=fileF)

        print("Matrix of core synteny scores [0,1] between genes in the family",file=fileF)
        printScoreMatrix(familyNum,subtreeL,familiesO,geneNames,scoresO,'coreSynSc',fileF)
        print(file=fileF)
        print(file=fileF)

        print("Matrix of synteny scores between genes in the family",file=fileF)
        printScoreMatrix(familyNum,subtreeL,familiesO,geneNames,scoresO,'synSc',fileF)
        print(file=fileF)
        print(file=fileF)

        printOutsideFamilyScores(familyNum,subtreeL,familiesO,geneNames,scoresO,fileF)
        print(file=fileF)
        print(file=fileF)

    def findIsland(searchStr,fileF=sys.stdout):
        '''Print the gene, family and island associated with searchStr. This
    is a wrapper that assumes various required objects are present at the
    top level.'''
        L=matchFamilyIsland(geneInfoD,geneNames,gene2FamD,fam2IslandD,searchStr)
        for geneName,fam,isl in L:
            print("<gene:"+str(geneName)+">","<family:"+str(fam)+">","<island:"+str(isl)+">",file=fileF)


    def printIsland(islandNum,synWSize,fileF=sys.stdout):
        '''Print the island and its genomic context in each species. We
        include synWSize/2 genes in either direction beyond the island.
        '''
        printIslandNeighb(islandNum,synWSize,subtreeL,islandByNodeL,familyL,geneOrderT,gene2FamD,fam2IslandD,geneInfoD,geneNames,strainNum2StrD,fileF)


    def printIslandsAtNode(nodeStr,fileF=sys.stdout):
        '''This is a wrapper to provide an easy way to print all the islands
    at a particular node in the tree. For ease of use, we take only a node
    number as argument, assuming all the other required stuff is available
    at the top level.
        '''
        node = strainStr2NumD[nodeStr]
        vPrintIslands(islandByNodeL[node],subtreeL,familyL,strainNum2StrD,geneNames,geneInfoD,fileF)


    def printCoreNonCoreByNode(fileF=sys.stdout):
        '''For each node in the focal clade, print the number of core and
    non-core families. That is, for that node, we get all families present
    in descendent species. We then look at their mrcas. Those with an mrca
    below the node in question are non-core, others are core.'''

        focTree = trees.subtree(tree,strainStr2NumD[paramD['rootFocalClade']])
        focInodesL=trees.iNodeList(focTree)
        familyByNodeL=createFamilyByNodeL(geneOrderT,gene2FamD)

        rowL=[]
        rowL.append(['Node','Core','Non-Core','Total','% Non-Core'])
        rowL.append(['----','----','--------','-----','----------'])
        for node in focInodesL:
            nonCore,core=coreNonCoreCtAtNode(tree,node,familyByNodeL,familyL)
            rowL.append([strainNum2StrD[node],str(core),str(nonCore),str(core+nonCore),str(format(nonCore/(core+nonCore),".3f"))])
        printTable(rowL,fileF=fileF)
        print(file=fileF)
        return


    ## Load data
    
    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)
    nodesL=trees.nodeList(tree)
    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])
    familiesO = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)
    # islandByNodeL=islands.readIslands(paramD['islandOutFN'],tree,strainStr2NumD)
    gene2FamD=createGene2FamD(familiesO)
    #fam2IslandD=createFam2IslandD(islandByNodeL)

    scoresO = scores.readScores(paramD['scoresFN'],geneNames)
    scoresO.createNodeConnectL(geneNames) # make nodeConnectL attribute
    # calc family error scores
    #families.calcErrorScores(familyL,scoresO,paramD['minNormThresh'],paramD['minCoreSynThresh'],paramD['minSynThresh'],paramD['famErrorScoreIncrementD'])
    
    code.interact(local=locals())

def debugWrapper(paramD):
    '''Take us into interactive mode for debugging.'''

    ## Set up the modules a bit differently for interactive mode
    import code,sys
    from .xenoGI import parameters,trees,genomes,families,islands,analysis,Score,scores

    import numpy
    from scipy.signal import find_peaks
    tree,strainStr2NumD,strainNum2StrD,geneNames,subtreeL,geneOrderT = loadMiscDataStructures(paramD)

    #scoresO=Score.Score()
    #scoresO.initializeDataAttributes(paramD['blastFilePath'],geneNames,strainStr2NumD)

    numBins = 80
    binWidth = 1.0/numBins # since scores range from 0-1
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)

    homologousPeakMissing,homologyRawThresholdD = scores.getHomologyRawThresholdD(scoresO)

    """
    scoreIterator = scoresO.iterateScoreByStrainPair((3,3),'rawSc')
    binHeightL,indexToBinCenterL = scores.scoreHist(scoreIterator,numBins)

    homologPeakWidth = 0
    widthRelHeight = 0.9
    homologRequiredProminence = 6
    homologLeftPeakLimit = 0.90
    homologRightPeakLimit = 1.0

    peakWidthInBins = homologPeakWidth / binWidth

    peakMinMax = [0, .05 / binWidth]

    tempBinHeightL = numpy.append(binHeightL,0)
    tempIndexToBinCenterL = numpy.append(indexToBinCenterL,1)

    peakIndL, propertiesD = find_peaks(tempBinHeightL, width = peakMinMax, rel_height = widthRelHeight, prominence = homologRequiredProminence)

    print(peakIndL)
    print(propertiesD)
    """
    
    code.interact(local=locals())

def simValidatonWrapper(paramD):
    '''Take us into interactive mode for running validation with
simulations. This will be a temporary flag only.'''

    ## Set up the modules a bit differently for interactive mode
    import code,sys
    from .xenoGI import parameters,trees,genomes,families,islands,analysis
    sys.path.insert(0,'/data/bushlab/htrans/dupDelProj/simCode')
    import validationSim

    moduleT=(parameters,trees,genomes,families,islands,analysis)
    
    xgiIslandGeneSetByNodeL,simIslandGeneSetByNodeL=validationSim.runValidation(moduleT,'simParams.py','params.py')
    
    code.interact(local=locals())
