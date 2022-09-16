import sys,numpy,os,random,glob,copy,shutil
sys.setrecursionlimit(100000)
from scipy.signal import find_peaks
from Bio import Phylo
from multiprocessing import Pool
from . import blast,genomes,trees,scores,DTLOR_DP,new_DTLOR_DP,islands
from .Family import *
from .Island import *
from .Tree import *
from .analysis import printTable

#### Globals

COSPECIATION = new_DTLOR_DP.NodeType.COSPECIATION
DUPLICATION =  new_DTLOR_DP.NodeType.DUPLICATION
LOSS =  new_DTLOR_DP.NodeType.LOSS
TRANSFER =  new_DTLOR_DP.NodeType.TRANSFER
LOCATION_ASSIGNMENT =  new_DTLOR_DP.NodeType.LOCATION_ASSIGNMENT
LOCATION_LIST =  new_DTLOR_DP.NodeType.LOCATION_LIST
LOCATION_MAPPING =  new_DTLOR_DP.NodeType.LOCATION_MAPPING
SPECIES_LIST =  new_DTLOR_DP.NodeType.SPECIES_LIST
SPECIES_MAPPING =  new_DTLOR_DP.NodeType.SPECIES_MAPPING
ORIGIN =  new_DTLOR_DP.NodeType.ORIGIN
ROOT =  new_DTLOR_DP.NodeType.ROOT
REARRANGEMENT =  new_DTLOR_DP.NodeType.REARRANGEMENT
ORIGIN_EVENT =  new_DTLOR_DP.NodeType.ORIGIN_EVENT

#### Main function

def createFamiliesO(speciesRtreeO,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,outputSummaryF):
    '''Main function to create families. First creates an initial families
object, then reconciles those against the species tree to make an
originFamilies object.
    '''

    # define some variables
    initFamilyFN = paramD['initFamilyFN']
    originFamilyFN =  paramD['originFamilyFN']
    geneInfoFN = paramD['geneInfoFN']
    numProcesses = paramD['numProcesses']
    D=int(paramD["duplicationCost"])
    T=int(paramD["transferCost"])
    L=int(paramD["lossCost"])
    O=int(paramD["originCost"])
    R=int(paramD["rearrangeCost"])

    # checks
    homologyCheck(genesO,aabrhHardCoreL,scoresO,outputSummaryF,paramD)

    # create blast families, output is directory of gene trees
    createBlastFamilies(paramD,speciesRtreeO,strainNamesT,scoresO,genesO,outputSummaryF)
    
    initialFamiliesO,locusMapD = createInitialFamiliesO(paramD,genesO,aabrhHardCoreL,scoresO,speciesRtreeO,outputSummaryF)
    
    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT)) # needed for writeFamilies
    writeFamilies(initialFamiliesO,initFamilyFN,genesO,strainNamesT,paramD)
    print("Initial families:",file=outputSummaryF)
    writeFamilyFormationSummary(initialFamiliesO,outputSummaryF)

    # reconcile
    initialFamiliesO = reconcileGeneTrees(initialFamiliesO.iterFamilies(),speciesRtreeO,initialFamiliesO,locusMapD,genesO,numProcesses,D,T,L,O,R)

    # reconcile for cases where family inserts repeatedly in same place
    if paramD['reconcilePermissiveOriginGeneListPath'] != None:
        initialFamiliesO = reconcilePermissiveOrigin(paramD,initialFamiliesO,speciesRtreeO,locusMapD,genesO)
    
    # create origin families
    initialFamiliesO,originFamiliesO = createOriginFamiliesO(speciesRtreeO,initialFamiliesO,paramD,genesO)
    print("Origin families:",file=outputSummaryF)
    writeFamilyFormationSummary(originFamiliesO,outputSummaryF)

    # write familes to file
    writeFamilies(initialFamiliesO,initFamilyFN,genesO,strainNamesT,paramD) # update with recons
    writeFamilies(originFamiliesO,originFamilyFN,genesO,strainNamesT,paramD)

    # optionally delete gene trees work dir
    if paramD['deleteGeneFamilyTreesDir']:
        shutil.rmtree(os.path.join(paramD['geneFamilyTreesDir']))
    
    return initialFamiliesO,originFamiliesO

## Support functions
    
def homologyCheck(genesO,aabrhHardCoreL,scoresO,outputSummaryF,paramD):
    '''Check if the number of genes in the hard core is low. Print warning if so. Also check for homology peak in raw scores histogram.'''

    ## Check number of hard core genes
    # figure out average number of genes per genome
    numGenes=0
    for rangeStart,rangeEnd in genesO.geneRangeByStrainD.values():
        if rangeEnd > numGenes:
            numGenes = rangeEnd
    avNumGenesPerGenome = (numGenes+1) / len(genesO.geneRangeByStrainD)

    propHC = round(len(aabrhHardCoreL) / avNumGenesPerGenome,3)
    
    if propHC < 0.2:
        print("Warning: the hard core gene set has only",len(aabrhHardCoreL),"genes. This is",propHC,file=outputSummaryF)
        print(""" of the average number of genes in the input genomes. A low number
 of core genes might result from one or more species being too
 distantly related, and could be associated with problems in family and
 island formation.""",file=outputSummaryF)

    ## Check for homology (right) peak in raw scores histogram
    scoreHistNumBins = paramD['scoreHistNumBins']
    binWidth = 1.0/scoreHistNumBins # since scores range from 0-1
    homologousPeakMissingL = []
    for strainPair in scoresO.getStrainPairs():

        if strainPair[0] == strainPair[1]:
            # strain vs self, skip
            continue
        
        scoreIterator = scoresO.iterateScoreByStrainPair(strainPair,'rawSc')
        binHeightL,indexToBinCenterL = scoreHist(scoreIterator,scoreHistNumBins)
        homologPeakLeftExtremePos=homologPeakChecker(binHeightL,indexToBinCenterL,binWidth,paramD)

        if homologPeakLeftExtremePos == float('inf'):
            homologousPeakMissingL.append(strainPair)
  
    if homologousPeakMissingL != []:
        print("""Warning: for one or more strain pairs, we failed to find a
 homologous (right) peak in the raw score histogram. The pair(s) in
 question are:""",file=outputSummaryF)
        for pairT in homologousPeakMissingL:
            print("  ",pairT[0]+'-'+pairT[1],file=outputSummaryF)
        
        print(""" A possible explanation for this failure is that one or more species
 is too distantly related. If this is the case it will result in poor
 family formation, and most families (e.g. 80% or more) will have genes
 in only one strain."""+"\n",file=outputSummaryF)

## Histograms and thresholds

def scoreHist(scoreIterator,scoreHistNumBins):
    '''Get a histogram with numpy, and return the bin height, and also a
list of indices to the middle position of each bin (in terms of the x value).'''
    binHeightL,edges = numpy.histogram(list(scoreIterator),bins=scoreHistNumBins,density=True)

    # make a list where the indices correspond to those of binHeightL,
    # and the values give the score value at the center of that bin
    indexToBinCenterL = []
    for i in range(1,len(edges)):
        center = (edges[i]+edges[i-1])/2
        indexToBinCenterL.append(center)

    return binHeightL,indexToBinCenterL

def homologPeakChecker(binHeightL,indexToBinCenterL,binWidth,paramD):
    '''Function to check for a peak due to homogology (right peak in
histogram). If such a peak exists, this function returns the position
(in score units) of the left most base of that peak. If no such peak
exits, this function returns infinity.

    '''
    peakL = [] # to collect in

    # in order to get a rightmost peak (if any) we add a dummy bin of
    # height 0 on right. add 1 to indexToBinCenterL for case of right
    # base of a peak in the last bin.
    tempBinHeightL = numpy.append(binHeightL,0)
    tempIndexToBinCenterL = numpy.append(indexToBinCenterL,1)
    
    # case 1 (normal case)
    L = findPeaksOneCase(tempBinHeightL,tempIndexToBinCenterL,binWidth,paramD['homologPeakWidthCase1'],paramD['widthRelHeight'],paramD['homologRequiredProminenceCase1'],paramD['homologLeftPeakLimitCase1'],paramD['homologRightPeakLimit'])
    peakL.extend(L)

    # case 2 (extreme prominence. But allow to be very narrow)
    L = findPeaksOneCase(tempBinHeightL,tempIndexToBinCenterL,binWidth,paramD['homologPeakWidthCase2'],paramD['widthRelHeight'],paramD['homologRequiredProminenceCase2'],paramD['homologLeftPeakLimitCase2'],paramD['homologRightPeakLimit'])
    peakL.extend(L)

    # case 3 (wide width with low prominence)
    L = findPeaksOneCase(tempBinHeightL,tempIndexToBinCenterL,binWidth,paramD['homologPeakWidthCase3'],paramD['widthRelHeight'],paramD['homologRequiredProminenceCase3'],paramD['homologLeftPeakLimitCase3'],paramD['homologRightPeakLimit'])
    peakL.extend(L)

    if peakL == []:
        return float('inf')
    else:
        # if there's more than one, we'll return the leftBasePos of the highest.
        peakL.sort(reverse=True)
        return peakL[0][2]

def findPeaksOneCase(binHeightL,indexToBinCenterL,binWidth,peakWidth,widthRelHeight,requiredProminence,leftPeakLimit,rightPeakLimit):
    '''Make one call to find_peaks. Peaks must be wider than peakWidth
(which is given in units of score) and more prominent than
requiredProminence, and fall between leftPeakLimit and
rightPeakLimit. Returns tuple
(peakHeight,peakPos,leftExtremeOfPeakPos,rightExtremeOfPeakPos) of any peaks
that meet criteria. All position values are returned in units of
score.

    '''
    peakWidthInBins = peakWidth / binWidth
    # we always measure width widthRelHeight down from peak toward base
    peakIndL, propertiesD = find_peaks(binHeightL, width = peakWidthInBins, rel_height = widthRelHeight, prominence = requiredProminence)

    # make sure they are to the right of leftPeakLimit
    peakPosInScoreUnitsL = []
    for i in range(len(peakIndL)):
        peakInd = peakIndL[i]
        peakHeight = binHeightL[peakInd]
        peakPos = indexToBinCenterL[peakInd]
        if leftPeakLimit < peakPos <= rightPeakLimit:
            # peak falls between the specified limits
            leftExtremeOfPeakPosInd = int(round(propertiesD["left_ips"][i]))
            leftExtremeOfPeakPos = indexToBinCenterL[leftExtremeOfPeakPosInd]
            rightExtremeOfPeakPosInd = int(round(propertiesD["right_ips"][i]))
            rightExtremeOfPeakPos = indexToBinCenterL[rightExtremeOfPeakPosInd]
            peakPosInScoreUnitsL.append((peakHeight,peakPos,leftExtremeOfPeakPos,rightExtremeOfPeakPos))
    return peakPosInScoreUnitsL

## Create blast families object

def createBlastFamilies(paramD,speciesRtreeO,strainNamesT,scoresO,genesO,outputSummaryF):
    '''Given a scoresO object, create gene families based on blast
    connectivity, then make gene trees with these.

    '''
    geneFamilyTreesDir = paramD['geneFamilyTreesDir']
    blastFamGeneTreeFileStem = paramD['blastFamGeneTreeFileStem']
    blastFamilyFN = paramD['blastFamilyFN']
    maxBlastFamSize = int(paramD['maxBlastFamSizeMultiplier'] * speciesRtreeO.leafCount())
    GeneRaxOutputDirN = paramD['GeneRaxOutputDirN']+"-"+blastFamGeneTreeFileStem
    
    ## get blast families as list of sets
    print("Blast families:",file=outputSummaryF)
    blastFamilySetL = createBlastFamilySetL(scoresO,genesO,strainNamesT,outputSummaryF,maxBlastFamSize)

    # save and add numbering so elements are (num,blastFamS)
    with open(blastFamilyFN,'w') as f:
        tempL=[] # for adding numbers
        i=1 # starting from 1 helps with some file name issues later
        for bSet in blastFamilySetL:
            f.write("\t".join(map(str,bSet))+"\n")
            tempL.append((i,bSet))
            i+=1
        blastFamilySetL = tempL
        del tempL
        
    ## make gene family trees for each blast family
    # create work dir if it doesn't already exist
    if glob.glob(geneFamilyTreesDir)==[]:
        os.mkdir(geneFamilyTreesDir)
        
    # delete any pre-existing blast family trees
    blastFamGeneTreeFilePath = os.path.join(geneFamilyTreesDir,blastFamGeneTreeFileStem+'*.tre')

    for fn in glob.glob(blastFamGeneTreeFilePath):
        os.remove(fn)

    # make trees
    if paramD['useGeneRaxToMakeSpeciesTrees']:
        # use GeneRax
        # run GeneRax on families of 3 genes or more. Run FastTree on rest
        # (will put branch lengths on 2 species trees, which we want
        geneRaxBlastFamilySetL = []
        fastTreeBlastFamilySetL = []
        for i,bSet in blastFamilySetL:
            if len(bSet) > 2:
                geneRaxBlastFamilySetL.append((i,bSet))
            else:
                fastTreeBlastFamilySetL.append((i,bSet))

        failedGeneRaxOrthoSetL = trees.makeGeneTreesGeneRax(paramD,False,genesO,geneFamilyTreesDir,blastFamGeneTreeFileStem,geneRaxBlastFamilySetL)

        fastTreeBlastFamilySetL.extend(failedGeneRaxOrthoSetL)
    else:
        # use FastTree for all
        fastTreeBlastFamilySetL = blastFamilySetL
        
    trees.makeGeneTreesFastTree(paramD,False,genesO,geneFamilyTreesDir,blastFamGeneTreeFileStem,fastTreeBlastFamilySetL)

    # collect and save blast family trees into a single file at the top level
    collectBlastFamilyTreesToSingleFile(paramD)

def createBlastFamilySetL(scoresO,genesO,strainNamesT,outputSummaryF,maxBlastFamSize):
    '''
    Input
    ------------------------------------------------
    scoresO:        a Score object that scores the 
                    BLAST results from all the genes

    Output
    ------------------------------------------------
    blastFamilySetL:      a list of sets where each set stores
                    all the genes that are connected as
                    indicated by significant BLAST score 
                    (that appear in scoresO)
    '''
   
    def getConnectedComponent(temp, gene, visited): 
        visited.add(gene)
        temp.add(gene) 
  
        # Repeat for all vertices adjacent 
        # to this gene
        neighbors=scoresO.getConnectionsGene(gene)
        if neighbors:
            for i in neighbors: 
                if i not in visited:     
                    # Update the list 
                    temp = getConnectedComponent(temp, i, visited) 
        
        return temp

    # main for createBlastFamilySetL
    if strainNamesT == None:
        allGenesL=list(genesO.iterGenes())
    else:
        allGenesL=list(genesO.iterGenes(strainNamesT))
    scoresO.createNodeConnectD() 
    connectedGenes=list(scoresO.nodeConnectD.keys())
    connecComponentSetL=[]
    visited=set()
    for gene in allGenesL: 
        if gene in scoresO.nodeConnectD:
            if gene not in visited: 
                temp =set()
                newFam=getConnectedComponent(temp, gene, visited)
                
                connecComponentSetL.append(newFam) 
        else:
            fam=set()
            fam.add(gene)
            connecComponentSetL.append(fam)

  
    # filter out largest
    blastFamilySetL,numAboveThresh = connecComponentSizeThreshold(connecComponentSetL,maxBlastFamSize,scoresO,outputSummaryF)

    summaryL = []
    summaryL.append(["Num initial connected components",str(len(connecComponentSetL))])
    summaryL.append(["Num components to split (have more than %d genes)"%maxBlastFamSize,str(numAboveThresh)])
    printTable(summaryL,indent=2,fileF=outputSummaryF)
    
    return blastFamilySetL

def connecComponentSizeThreshold(connecComponentSetL,maxBlastFamSize,scoresO,outputSummaryF):
    '''For family sets larger than maxBlastFamSize, split by raising the
threshold for homology and reclustering.'''

    blastFamilySetL = []
    tooBigL = []
    numAboveThresh = 0
    for clusterS in connecComponentSetL:
        if len(clusterS)> maxBlastFamSize:
            numAboveThresh+=1
            # force split
            tooBigL.append((len(clusterS),clusterS))
        else:
            blastFamilySetL.append((len(clusterS),clusterS))
            
    # split the big ones
    tooBigL.sort(reverse=True,key=lambda x: x[0])
    for _,fullClusterS in tooBigL:
        
        splitClusterL = splitClusterFailsafe(fullClusterS,scoresO,maxBlastFamSize)

        for clusterS in splitClusterL:
            blastFamilySetL.append((len(clusterS),clusterS))

    blastFamilySetL.sort(reverse=True,key=lambda x: x[0])

    # strip out length and return
    return [clusterS for _,clusterS in blastFamilySetL],numAboveThresh

def splitClusterFailsafe(fullClusterS,scoresO,maxBlastFamSize):
    '''Split fullClusterS by setting higher and higher thresholds for
homology. Find a threshold that is as low as possible, but still
splits fullClusterS such that no subclusters are larger than
maxBlastFamSize.

    '''
    splitClusterBinSearchSteps = 25 # max value of steps we'll allow
    
    connecT = getConnectionT(fullClusterS,scoresO)

    lind = 0
    hind = len(connecT)-1

    bestInd = splitBinSearch(connecT,lind,hind,fullClusterS,maxBlastFamSize,splitClusterBinSearchSteps)
    splitClusterL = buildClusterLFromConnecT(connecT[bestInd:],fullClusterS)
   
    return splitClusterL
    
def getConnectionT(fullClusterS,scoresO):
    '''Get a tuple representation of the network for this cluster. Each
    element is ((node1,node2),rawSc), and we sort of rawSc.
    '''
    connecD={}
    for gene1 in fullClusterS:
        neighborL=scoresO.getConnectionsGene(gene1)
        for gene2 in neighborL:
            if gene1 != gene2:
                # must eliminate self connections
                keyT=tuple(sorted([gene1,gene2]))
                connecD[keyT] = scoresO.getScoreByEndNodes(keyT[0],keyT[1],'rawSc')
    connecT = tuple(sorted(connecD.items(),key=lambda x: x[1]))
    return connecT

def splitBinSearch(connecT,lind,hind,fullClusterS,maxBlastFamSize,step):
    '''Find the lowest index into connecT where the network made from that
point forward has no clusters larger than maxBlastFamSize.

    '''
    mind = int((lind+hind)/2)
    # base case
    if mind==lind or mind==hind:
        return hind
    elif step == 0:
        return hind
    else:
        splitClusterL = buildClusterLFromConnecT(connecT[mind:],fullClusterS)
        if isSplitClusterLSmallEnough(splitClusterL,maxBlastFamSize):
            hind=mind
        else:
            lind=mind
        # recurse.
        return splitBinSearch(connecT,lind,hind,fullClusterS,maxBlastFamSize,step-1)
        
def isSplitClusterLSmallEnough(splitClusterL,maxBlastFamSize):
    '''Test if the clusters in splitClusterL are below
maxBlastFamSize. Return True if they all are, False if any are too
big.
    '''
    for clusterS in splitClusterL:
        if len(clusterS)>maxBlastFamSize:
            return False
    return True
    
def buildClusterLFromConnecT(connecT,fullClusterS):
    '''Build the clusters implied by connecL (which has already been
reduced by thresholding). fullClusterS gives the full set of genes we
must divide.

    '''
    ## functions
    def getNodeConnectD(connecT):
        '''Take all the connections in connecT, and put them in a node based
    representaiton in nodeConnectD.'''

        nodeConnectD = {}
        for (gn1,gn2),_ in connecT:

            if gn1 not in nodeConnectD:
                nodeConnectD[gn1] = set([gn2])
            else:
                nodeConnectD[gn1].add(gn2)

            if gn2 not in nodeConnectD:
                nodeConnectD[gn2] = set([gn1])
            else:
                nodeConnectD[gn2].add(gn1)

        return nodeConnectD

    def getConnectedComponent(visitedS,clusterS,nodeConnectD,gene):
        '''Get all the genes connected to gene in nodeConnectD'''
        visitedS.add(gene)
        clusterS.add(gene) 

        # Repeat for all vertices adjacent to this gene
        neighborsL= nodeConnectD[gene]
        for neighbGene in neighborsL: 
            if neighbGene not in visitedS:     
                # Update the list 
                clusterS = getConnectedComponent(visitedS,clusterS,nodeConnectD,neighbGene) 

        return clusterS

    ## buildClusterLFromConnecT main
    nodeConnectD = getNodeConnectD(connecT)
    
    splitClusterL = [] # for output
    visitedS=set()
    for gene in fullClusterS: 
        if gene in nodeConnectD:
            if gene not in visitedS: 
                tempClusterS=set()
                tempClusterS=getConnectedComponent(visitedS,tempClusterS,nodeConnectD,gene)
                
                splitClusterL.append(tempClusterS) 
        else:
            fam=set()
            fam.add(gene)
            splitClusterL.append(fam)
        
    # check if there is redundancy
    L=[]
    S=set()
    for clusterS in splitClusterL:
        L.extend(clusterS)
        S.update(clusterS)

    if len(L) != len(S):
        print("inner---")
        print(len(L))
        print(len(S))

    return splitClusterL


def addGeneToSplitClusterL(splitClusterL,gn1,gn2):
    '''Given a connection pair, gn1,gn2, add to connecL, creating a new
set in connecL if none of the existing sets have either gene.

    '''
    for clusterS in splitClusterL:
        if gn1 in clusterS:
            clusterS.add(gn2)
            return splitClusterL
        if gn2 in clusterS:
            clusterS.add(gn1)
            return splitClusterL
    # if we get here, must create new cluster
    newClusterS = set([gn1,gn2])
    splitClusterL.append(newClusterS)
    return splitClusterL

def collectBlastFamilyTreesToSingleFile(paramD):
    '''Collect and save blast family trees into a single file at the top
level.
    '''
    blastFamGeneTreeL = loadGeneTreesFromDir(paramD,paramD['blastFamGeneTreeFileStem'])
    
    with open(paramD['blastFamilyTreeFN'],'w') as f:
        for famNum,geneUtreeO in blastFamGeneTreeL:
            treeStr = geneUtreeO.toNewickStr(includeBrLength=True)
            f.write(str(famNum)+'\t'+treeStr+'\n')
    return

## Initial family formation    
def createInitialFamiliesO(paramD,genesO,aabrhHardCoreL,scoresO,speciesRtreeO,outputSummaryF):
    ''''''
    blastFamGeneTreeFileStem = paramD['blastFamGeneTreeFileStem']
    aabrhHardCoreGeneTreeFileStem = paramD['aabrhHardCoreGeneTreeFileStem']
    geneFamilyTreesDir = paramD['geneFamilyTreesDir']
    maxInitialFamSize = int(paramD['maxInitialFamSizeMultiplier'] * speciesRtreeO.leafCount())
    forceSplitUtreeBalanceMultiplier = paramD['forceSplitUtreeBalanceMultiplier']

    ## load gene trees
    blastFamGeneTreeL = loadGeneTreesFromDir(paramD,blastFamGeneTreeFileStem)

    # delete any aabrh gene trees
    # aside: for aabrh trees made for species tree calculation, the
    # headers are wrong, containing the species rather than the
    # gene. These should be in a different directory anyway...


    aabrhGtFilePath = os.path.join(geneFamilyTreesDir,aabrhHardCoreGeneTreeFileStem+'*.tre')
    aabrhGtFileNameL=list(glob.glob(aabrhGtFilePath))
    for fn in aabrhGtFileNameL:
        os.remove(fn)
    newAabrhHardCoreL = [] # add numbering
    orthoNum = 1 # start at 1
    for orthoT in aabrhHardCoreL:
        newAabrhHardCoreL.append((orthoNum,orthoT))
        orthoNum += 1

    trees.makeGeneTreesFastTree(paramD,False,genesO,geneFamilyTreesDir,aabrhHardCoreGeneTreeFileStem,newAabrhHardCoreL)
    
    aabrhHardCoreGeneTreeL = loadGeneTreesFromDir(paramD,aabrhHardCoreGeneTreeFileStem)

    ## split blast family trees

    # see how many in blastFamGeneTreeL match aabrhHardCoreGeneSet
    # get sets of aabrh genes so we can avoid splits in these families
    aabrhHardCoreGeneSetL = [set(aabrhUtreO.leaves()) for famNum,aabrhUtreO in aabrhHardCoreGeneTreeL]
    
    # obtain threshold. Will split on branches larger than splitThresh
    splitThresh = calculateTreeSplitThreshold(paramD,aabrhHardCoreGeneTreeL)

    singleGeneTreeL = [] # need this since splitting funcs don't like single tip trees
    geneTreeL = []
    for blastFamNum,blastFamUtreeO in blastFamGeneTreeL:
        if blastFamUtreeO.leafCount() == 1:
            singleGeneTreeL.append((blastFamUtreeO.leafCount(),blastFamNum,blastFamUtreeO))
        else:

            # pull out aabrh sets that are subsets of this tree. We
            # will split so as to ensure that the genes of the aabrh
            # set are not separated.
            subsetAabrhL = getSubsetAabrhL(set(blastFamUtreeO.leaves()),aabrhHardCoreGeneSetL)
            
            for aUtO in splitUtreeThreshold([blastFamUtreeO],splitThresh,subsetAabrhL,set()):
                if aUtO.leafCount() == 1:
                    singleGeneTreeL.append((aUtO.leafCount(),blastFamNum,aUtO))
                else:
                    # failsafe to keep below maxInitialFamSize. If none too
                    # big, will pass back same list

                    for bUtO in splitUtreeFailsafe([aUtO],maxInitialFamSize,forceSplitUtreeBalanceMultiplier,subsetAabrhL):
                        geneTreeL.append((bUtO.leafCount(),blastFamNum,bUtO))

    geneTreeL.extend(singleGeneTreeL)
    del singleGeneTreeL
    geneTreeL.sort(reverse=True,key=lambda x: x[0]) # big families first
    
    ## make initial families object
    
    # create output objects
    locusMapD={}
    # construct Family object for initial families
    initialFamiliesO = Families(speciesRtreeO)

    # get synteny thesholds for locus family formation
    synThresholdD = getSynThresholdD(paramD,scoresO,genesO,aabrhHardCoreL,speciesRtreeO)

    # store in familiesO object
    # these locusFam Ids need to be unique across all, start counting from 1.
    famNumCounter=1
    locusFamNumCounter=1
    totalAddedTolocusFamilies=0

    for size,blastFamNum,utO in geneTreeL:
        familyS = set(map(int,utO.leaves())) # convert from string to numeric
        locusFamLL=divideInitialFamilyIntoLocusFamilies(familyS,genesO,scoresO,paramD,synThresholdD)
        locusFamNumCounter,famNumCounter,totalAddedTolocusFamilies = addFamilyToInitialFamiliesO(familyS,initialFamiliesO,famNumCounter,locusFamLL,locusMapD,locusFamNumCounter,genesO,speciesRtreeO,totalAddedTolocusFamilies,utO,blastFamNum)
        
    return initialFamiliesO, locusMapD

def loadGeneTreesFromDir(paramD,geneTreeFileStem):
    '''Load gene trees from the geneFamilyTreesDir that begin with
    geneTreeFileStem.
    '''
    geneFamilyTreesDir = paramD['geneFamilyTreesDir']
    if not os.path.isdir(geneFamilyTreesDir):
        raise FileNotFoundError("Directory of gene trees is missing.")

    # load gene trees, divide into bifurcating vs. multifurcating
    allGtFilePath = os.path.join(geneFamilyTreesDir,geneTreeFileStem+'*.tre')
    allTreeFN_L=list(sorted(glob.glob(allGtFilePath)))

    geneTreeL = []
    for geneTreeFN in allTreeFN_L:
        geneUtreeO = Utree()
        geneUtreeO.fromNewickFile(geneTreeFN)
        famNum = int(geneTreeFN.split(geneTreeFileStem)[1].split('.tre')[0].lstrip('0'))
        geneTreeL.append((famNum,geneUtreeO))

    return geneTreeL

def getSubsetAabrhL(famS,aabrhHardCoreGeneSetL):
    '''Identify aabrh sets that are subsets of famS. Return list of sets.
    '''
    subsetAabrhL = []
    for aabrhS in aabrhHardCoreGeneSetL:
        if aabrhS.issubset(famS):
            subsetAabrhL.append(aabrhS)

    return subsetAabrhL

def calculateTreeSplitThreshold(paramD,aabrhHardCoreGeneTreeL):
    '''Given a list of aabrh hard core gene trees, calculate a threshold
for splitting gene trees. For each aabrh tree, we get the maximum
branch lenth. Then based on this distribution of lengths we obtain a
threshold.

    '''
    multiplierForObtainingSplitThresholds = paramD['multiplierForObtainingSplitThresholds']
    quantileForObtainingSplitThresholds = paramD['quantileForObtainingSplitThresholds']
    
    maxBranchLenL = []
    for _,geneUtreeO in aabrhHardCoreGeneTreeL:
        _,brLen = geneUtreeO.maxBranchLen()
        maxBranchLenL.append(brLen)
    
    splitThresh = multiplierForObtainingSplitThresholds * numpy.quantile(maxBranchLenL,quantileForObtainingSplitThresholds)

    return splitThresh

def splitUtreeThreshold(utreeL,splitThresh,subsetAabrhL,doNotSplitBranchPairS):
    '''Given an input list containing Utree objects, recursively split
these until there are no branches longer than splitThresh. Trees in
utreeL should not be single tip trees. 
    '''

    def getBranchToSplit(utreeL,doNotSplitBranchPairS):
        i=-1 # for empty list case
        for i in range(len(utreeL)):
            for splitBranchPair,splitBrLen in utreeL[i].getBranchesByLengthL():
                if splitBrLen > splitThresh and splitBranchPair not in doNotSplitBranchPairS:
                    return i,splitBranchPair,splitBrLen
                elif splitBrLen < splitThresh:
                    # advance to next tree, all remaining branches in this one below thresh
                    break
        return i+1,None,None # none above thresh
        
    ## main splitUtree

    # get tree index in utreeL and branch to split
    aboveThreshInd,splitBranchPair,splitBrLen = getBranchToSplit(utreeL,doNotSplitBranchPairS)

    if splitBranchPair == None: # done
        return utreeL
    else:
        completeL = utreeL[:aboveThreshInd] # these are done

        # split the branch we found that was too large. aboveThreshInd
        # was chosen to avoid stuff in doNotSplitBranchPairS

        remainderL = []
        aUtreeO,bUtreeO = utreeL[aboveThreshInd].split(splitBranchPair)
        
        # check that this split preserves aabrh sets
        if dividesAabrhSets(aUtreeO,bUtreeO,subsetAabrhL):
            doNotSplitBranchPairS.add(splitBranchPair)
            remainderL.append(utreeL[aboveThreshInd])
        else:
            # did not divide aabrh sets. Let's go with it.
            # if these split products have one tip, put in completeL,
            # otherwise in remainderL
            if aUtreeO.leafCount() == 1:
                completeL.append(aUtreeO)
            else:
                remainderL.append(aUtreeO)

            if bUtreeO.leafCount() == 1:
                completeL.append(bUtreeO)
            else:
                remainderL.append(bUtreeO)

            # When we split off utrees, we merge 2 branches into a
            # larger branch which should not be split unless one of
            # the original branches that went into it is above
            # theshold
            doNotSplitBranchPairS = checkIfNewlyMergedBranchShouldNotBeSplit(utreeL[aboveThreshInd],splitThresh,splitBranchPair,doNotSplitBranchPairS)
                
        # put the results back in the todo list and recurse
        remainderL.extend(utreeL[aboveThreshInd+1:])
        remainderSplitL = splitUtreeThreshold(remainderL,splitThresh,subsetAabrhL,doNotSplitBranchPairS)
        completeL.extend(remainderSplitL)
        return completeL

def dividesAabrhSets(aUtreeO,bUtreeO,subsetAabrhL):
    '''Returns True if one or more of the sets in subsetAabrhL are split
accross aUtreeO and bUtreeO. Returns False otherwise. In this case the
sets are a subset of aUtreeO, a subset of bUtreeO or a subset of
neither.
    '''
    aS = set(aUtreeO.leaves())
    bS = set(bUtreeO.leaves())
    for aabrhS in subsetAabrhL:
        if aabrhS.intersection(aS) and aabrhS.intersection(bS):
            # overlap with both. aabrh has been divided.
            return True
    return False

def checkIfNewlyMergedBranchShouldNotBeSplit(origUtree,splitThresh,splitBranchPair,doNotSplitBranchPairS):
    '''In making utrees on either side of the split, we have removed two
nodes. Where a node is removed, two branches are merged. this bigger
merged branch should only be split later if one or both of the branches
which went into it were above threshold.
    '''

    # funcs

    def somethingAboveThreshold(sbtree,node,splitThresh):
        '''Determine if any of the children of the 2 subtrees are above threshold.'''
        # get subtrees on each side
        for child in sbtree.children(node):
            if sbtree.getBranchLen((node,child)) > splitThresh:
                return True
        return False

    
    # main section of checkIfNewlyMergedBranchShouldNotBeSplit
    # root orig tree on splitBranchPair
    rtree = origUtree.rootIncludeBranchLen(splitBranchPair)

    for node in splitBranchPair:
        sbtree = rtree.subtree(node)
    
        if not somethingAboveThreshold(sbtree,node,splitThresh):
            # the children of this subtree define a branch in the
            # split utree. We should not split on this branch, since
            # neither of the branches in the rooted subtree are above
            # threshold.
            branchNotToSplitPair = sbtree.children(node)
            if len(branchNotToSplitPair) == 2:
                
                doNotSplitBranchPairS.add(branchNotToSplitPair)
                doNotSplitBranchPairS.add((branchNotToSplitPair[1],branchNotToSplitPair[0])) # both ways
            
    return doNotSplitBranchPairS
    
    
def splitUtreeFailsafe(utreeL,maxInitialFamSize,forceSplitUtreeBalanceMultiplier,subsetAabrhL):
    '''Function for cutting tree size down in order to limit dtlor
calculation time.  Given an input list containing Utree objects,
recursively split these until all are below maxInitialFamSize. Trees in
utreeL should not be single tip trees.
    '''

    def getIndFirstTooBig(utreeL):
        i=-1 # for empty list case
        for i in range(len(utreeL)):
            if utreeL[i].leafCount() > maxInitialFamSize:
                return i
        return i+1 # none too big
        
    # main splitUtreeFailsafe
    tooBigInd = getIndFirstTooBig(utreeL)

    if tooBigInd == len(utreeL): # done
        return utreeL
    else:
        completeL = utreeL[:tooBigInd] # these are done
        
        # split the tree that was too large
        aUtreeO,bUtreeO = forceSplitUtree(utreeL[tooBigInd],forceSplitUtreeBalanceMultiplier,subsetAabrhL)
        
        remainderL = []

        # if these split products have one tip, put in completeL,
        # otherwise in remainderL
        if aUtreeO.leafCount() == 1:
            completeL.append(aUtreeO)
        else:
            remainderL.append(aUtreeO)

        if bUtreeO.leafCount() == 1:
            completeL.append(bUtreeO)
        else:
            remainderL.append(bUtreeO)
            
        # put the rest of utreeL in the to do list and recurse
        remainderL.extend(utreeL[tooBigInd+1:])
        remainderSplitL = splitUtreeFailsafe(remainderL,maxInitialFamSize,forceSplitUtreeBalanceMultiplier,subsetAabrhL)
        completeL.extend(remainderSplitL)
        return completeL

def forceSplitUtree(geneUtreeO,forceSplitUtreeBalanceMultiplier,subsetAabrhL):
    '''Helper function for splitting excessively large trees (in order to
cut dtlor calculation time). Given a large unrooted gene tree, split
on an internal branch and return the two resulting subtrees. We use
two criteria to pick branches. 1. branch length (larger is better)
2. the balance of leaves on the left and right sides. (more balanced
is better). We integrate these by sorting the list of branches on each
metric independently. We then minimize the sum of the indices in the
two lists.

    '''
    # get branches sorted by branchLen, high to low
    branchLenL = geneUtreeO.getBranchesByLengthL()
    branchIndD = {}
    i=0
    for branchPair,brLen in branchLenL:
        branchIndD[branchPair] = i
        i+=1

    # branches sorted by node balance
    balanceL = []
    for branchPair,brLen in branchLenL:
        bal = branchBalanceCalc(branchPair,geneUtreeO)
        balanceL.append((branchPair,bal))
    balanceL.sort(key=lambda x: x[1])
    branchBalanceD = {}
    i=0
    for branchPair,bal in balanceL:
        branchBalanceD[branchPair] = i
        i+=1

    # combine
    comboL = []
    for branchPair in geneUtreeO.branchPairT:
        comboVal = branchIndD[branchPair] + forceSplitUtreeBalanceMultiplier * branchBalanceD[branchPair]
        comboL.append((branchPair,comboVal))
    comboL.sort(key=lambda x: x[1]) # lower comboVals better

    # split. We assume there is always some branch where we can split and not divide aabrh's.
    for splitBranchPair,_ in comboL:
        aUtreeO,bUtreeO = geneUtreeO.split(splitBranchPair)
        if not dividesAabrhSets(aUtreeO,bUtreeO,subsetAabrhL):
            # these two trees do not divide aabrh's
            break

    return aUtreeO,bUtreeO

def branchBalanceCalc(branchPair,geneUtreeO):
    '''Calculate the number of leaves to the left of branchPair (l) and to
the right of branchPair (r). return abs(l-r). The smaller this value,
the more we want to split on a branch.
    '''
    aO,bO = geneUtreeO.split(branchPair)
    l=aO.leafCount()
    r=bO.leafCount()
    return abs(l-r)

def getSynThresholdD(paramD,scoresO,genesO,aabrhHardCoreL,speciesRtreeO):
    '''Creates a dictionary to store synteny thresholds. This dictionary
itself contains two dictionaries one each for minCoreSynThresh
(minimum core synteny score allowed for family formation), and
minSynThresh (minimum synteny score allowed for family
formation). These dictionaries in turn are keyed by strain pair.
    '''
    quantileForObtainingSynThresholds = paramD['quantileForObtainingSynThresholds']
    multiplierForObtainingSynThresholds = paramD['multiplierForObtainingSynThresholds']
    
    synThresholdD = {}
    synThresholdD['minSynThreshold'] = {}
    synThresholdD['minCoreSynThreshold'] = {}
    
    # coreSynSc. Note that if this ends up higher than 0.5, then there
    # will need to be at least one core gene on both sides.
    for strainPair in scoresO.getStrainPairs():
        aabrhScL = scores.getScoresStrainPair(scoresO,strainPair,'coreSynSc',genesO,aabrhHardCoreL)
        thresh = multiplierForObtainingSynThresholds * numpy.quantile(aabrhScL,quantileForObtainingSynThresholds)
        synThresholdD['minCoreSynThreshold'][strainPair] = thresh

    # synSc
    for strainPair in scoresO.getStrainPairs():
        aabrhScL = scores.getScoresStrainPair(scoresO,strainPair,'synSc',genesO,aabrhHardCoreL)
        thresh = multiplierForObtainingSynThresholds * numpy.quantile(aabrhScL,quantileForObtainingSynThresholds)
        synThresholdD['minSynThreshold'][strainPair] = thresh
        
    # In the case of family formation at a tip, we're interested in
    # genes that duplicated after the last species split off. So the
    # thresholds at a tip really shouldn't be based on the given
    # genome against itself. Basing them instead on what we saw at the
    # parent node (of the last divergence) seems reasonable, and is
    # what we'll do here. We'll now replace the entries in
    # synThresholdD for tips with the values from the parent node.

    for strainPair in scoresO.getStrainPairs():
        if strainPair[0] == strainPair[1]:
            # Tip
            leafStrain = strainPair[0]
            synThresholdD['minCoreSynThreshold'][strainPair] = getTipThreshold(speciesRtreeO,leafStrain,synThresholdD,'minCoreSynThreshold')
            synThresholdD['minSynThreshold'][strainPair] = getTipThreshold(speciesRtreeO,leafStrain,synThresholdD,'minSynThreshold')

    return synThresholdD

def getTipThreshold(speciesRtreeO,leafStrain,synThresholdD,thresholdType):
    '''Get the synteny threshold values at a tip called leafStrain by
averaging the values between that strain and its nearest neighbors.'''
    
    neighbL = speciesRtreeO.getNearestNeighborL(leafStrain)
    avThresh = 0
    for neighb in neighbL:
        strainPair = tuple(sorted((neighb,leafStrain)))
        avThresh += synThresholdD[thresholdType][strainPair]
    avThresh /= len(neighbL)
    return avThresh

def addFamilyToInitialFamiliesO(familyS,initialFamiliesO,famNumCounter,locusFamLL,locusMapD,locusFamNumCounter,genesO,speciesRtreeO,totalAddedTolocusFamilies,geneUtreeO,sourceFam):
    '''Given a families object, a family set and a list of locus families
in that family, add the family to the object.
    '''
    speciesL=list(set([genesO.numToStrainName(gene) for gene in familyS]))
    mrca=speciesRtreeO.findMrca(speciesL)
    # add each initial family as a Family object (still empty)
    initialFamiliesO.initializeFamily(famNumCounter,mrca,"initial",geneTreeO=geneUtreeO,sourceFam=sourceFam)
    for locusFamilyL in locusFamLL: #locusFamilyL contains all the genes in that lf
        speciesL=[]
        for gene in locusFamilyL:
            locusMapD[gene]=locusFamNumCounter
            speciesL.append(genesO.numToStrainName(gene))
        lfMrca=speciesRtreeO.findMrca(speciesL)
        # create a lf. for initial families, locusNum is the same as locusFamNum.
        lf=LocusFamily(famNumCounter,locusFamNumCounter,lfMrca,locusNum=locusFamNumCounter)
        lf.addGenes(locusFamilyL, genesO)
        totalAddedTolocusFamilies+=len(locusFamilyL)
        initialFamiliesO.addLocusFamily(lf)
        locusFamNumCounter+=1
    famNumCounter+=1
    return locusFamNumCounter,famNumCounter,totalAddedTolocusFamilies

def divideInitialFamilyIntoLocusFamilies(familyS, genesO, scoresO, paramD,synThresholdD):
    """
    Input 
    ----------------------------------------
    familyS:            An initial family of genes generated by 'createBlastFamilySetL'. This is represented 
                        by a set that contains genes that are connected by BLAST results
    scoresO:            The Score object that contains the information for all pairwise synteny scores

    Return 
    ----------------------------------------
    locusFamilies:      A list of lists, sub-lists contain locus families produced by clustering

    """

    minGenes = 2 # don't divide if num genes less than this.

    
    #if there are too few genes, we should probably not cluster further
    if len(familyS) < minGenes:
        locusFamilyL=list(familyS)
        return [locusFamilyL]

    # construct the fully connected graph with weighted edges 
    num_gene=len(familyS)
    genesL=list(familyS)
    geneNeighborsD={}
    # there will be no edge going from one to self.

    # initialize
    for gene in genesL:
        geneNeighborsD[gene] = []
    
    #loop over all pairs of genes
    for i in range(len(genesL)):
        for j in range(i+1, len(genesL)):
            gene1=genesL[i]
            gene2=genesL[j]
            if scoresO.isEdgePresentByEndNodes(gene1,gene2):
                sameLocusFam = isSameLocusFamily(gene1,gene2,scoresO,genesO,paramD,synThresholdD)
                # the graph is symmetrical

                # if above threshold for both syn sc, add as neighbors
                if sameLocusFam:
                    geneNeighborsD[gene1].append(gene2)
                    geneNeighborsD[gene2].append(gene1)

    # find connected components
    def getConnectedComponent(tempFamL, gene, visitedS,geneNeighborsD):
        '''Find components connected to gene and put in tempFamL.'''
        visitedS.add(gene)
        tempFamL.append(gene)
        neighborsL=geneNeighborsD[gene]
        if len(neighborsL)>0:
            for neighbGene in neighborsL: 
                if neighbGene not in visitedS:     
                    # Update the list 
                    tempFamL = getConnectedComponent(tempFamL, neighbGene, visitedS, geneNeighborsD) 
        return tempFamL 

    locusFamLL=[]
    visitedS=set()
    for gene in genesL: 
        if gene not in visitedS: 
            tempFamL =[]
            newFamL=getConnectedComponent(tempFamL, gene, visitedS, geneNeighborsD)
            locusFamLL.append(newFamL) 

    return locusFamLL

def isSameLocusFamily(gene1,gene2,scoresO,genesO,paramD,synThresholdD):
    '''Given two genes in the same family, determine if they meet the
synteny requirements to be put in the same LocusFamily. Returns
boolean.
    '''
    if not scoresO.isEdgePresentByEndNodes(gene1,gene2):
        # Within our families, there may be some gene-gene edges
        # missing due to the fact that blast could have just missed
        # significance etc. If the edge isn't there, then we do not
        # have evidence that these genes should be in the same locus
        # family, and we return false.
        
        # NOTE: An alternative approach would be to actually calculate
        # these scores here. But they're likely to be low...

        return False
    
    coreSynSc = scoresO.getScoreByEndNodes(gene1,gene2,'coreSynSc')
    synSc = scoresO.getScoreByEndNodes(gene1,gene2,'synSc')

    strain1 = genesO.numToStrainName(gene1)
    strain2 = genesO.numToStrainName(gene2)
    strainPair = tuple(sorted([strain1,strain2]))
    minSynThreshold = synThresholdD['minSynThreshold'][strainPair]
    minCoreSynThreshold = synThresholdD['minCoreSynThreshold'][strainPair]
    
    if coreSynSc < minCoreSynThreshold or synSc < minSynThreshold:
        # one of the two types of synteny below threshold, so this
        # pair doesn't meet the requirements for being in the same
        # LocusFamily
        addIt = False
    else:
        addIt = True
            
    return addIt

#### Family formation summary 

def writeFamilyFormationSummary(familiesO,outputSummaryF):
    '''Write a summary of the families in familiesO to the file handle
outputSummaryF.
    '''

    summaryL=[]
    summaryL.append(["Total number of Families",str(len(familiesO.familiesD))])
    summaryL.append(["Total number of LocusFamilies",str(len(familiesO.locusFamiliesD))])
    summaryL.append(["Total number of genes in families",str(len(familiesO.getAllGenes()))])

    singleLfFams = 0
    multipleLfFams = 0
    singleStrainFams = 0
    multipleStrainFams = 0
    famLenL = []
    for fam in familiesO.iterFamilies():
        famLenL.append(fam.geneCount())
        # locus families
        if len(fam.getLocusFamilies()) == 1:
            singleLfFams += 1
        else:
            multipleLfFams += 1
        # strains
        if len(set(fam.iterStrains())) == 1:
            singleStrainFams += 1
        else:
            multipleStrainFams += 1

    sizeDistSummaryL = [numpy.quantile(famLenL,0),numpy.quantile(famLenL,.25),numpy.quantile(famLenL,.5),numpy.quantile(famLenL,.75),numpy.quantile(famLenL,1)]
    summaryL.append(["Family size dist [min,q1,q2,q3,max]",str(sizeDistSummaryL)])
    
    summaryL.append(["Number of families with one LocusFamily",str(singleLfFams)])
    summaryL.append(["Number of families with multiple LocusFamilies",str(multipleLfFams)])
    summaryL.append(["Number of families with gene(s) in only one strain",str(singleStrainFams)])
    summaryL.append(["Number of families with genes in multiple strains",str(multipleStrainFams)])
    
    printTable(summaryL,indent=2,fileF=outputSummaryF)

    return

#### Reconciliation

def reconcileGeneTrees(ifamIter,speciesRtreeO,initialFamiliesO,locusMapD,genesO,numProcesses,D,T,L,O,R):
    '''Reconcile gene family trees to the species tree using the DTLOR
algorithm. ifamIter is an iterator (or list) of initial families to
work on. D,T,L,O,R are the DTLOR costs.
    '''

    argumentL = []
    for iFamO in ifamIter:
        initFamNum = iFamO.famNum
        geneUtreeO = iFamO.geneTreeO

        # skip single gene
        if iFamO.geneCount() == 1:
            continue

        # in loop
        tipMapD=getTipMapping(geneUtreeO,genesO)

        # Make new gtLocusMapD with only those genes in geneUtreeO
        gtLocusMapD = reduceLocusMap(geneUtreeO,locusMapD)

        # ensure tree is binary
        if len(geneUtreeO.multifurcatingNodes()) > 0:
            # it's multifurcating, arbitrarily binarize
            # In future, we can implement dtlor for multifurcating nodes
            geneUtreeO = geneUtreeO.binarize(gtLocusMapD)
      
        # add to argumentL
        argT = (initFamNum,speciesRtreeO,geneUtreeO,tipMapD,gtLocusMapD,D,T,L,O,R)
        argumentL.append(argT)

    # run on multiple processors
    with Pool(processes=numProcesses) as p:
        for initFamNum,optGeneRtreeO,optG,minCost in p.imap_unordered(reconcileOneUnRootedGeneTree, argumentL):
            
            # store
            ifam = initialFamiliesO.getFamily(initFamNum)
            ifam.addGeneTree(optGeneRtreeO)
            ifam.addGraphD(optG)
            ifam.dtlorCost = minCost

    return initialFamiliesO
    
def getTipMapping(geneUtreeO, genesO):
    """
    Fill out the tip mapping (from gene to species) using the function from genomes.
    """
    tipMapD={}
    for leaf in geneUtreeO.leaves():
        # the leaf is a gene number in string form
        tipMapD[leaf]=genesO.numToStrainName(int(leaf))
    return tipMapD

def reduceLocusMap(geneUtreeO,locusMapD):
    '''Create a new locus map D with only entries for genes in geneUtreeO.'''
    gtLocusMapD={}
    for leaf in geneUtreeO.leaves():
        # the leaf is a gene number in string form
        gtLocusMapD[leaf] = locusMapD[int(leaf)]
    return gtLocusMapD
        
def reconcileOneUnRootedGeneTree(argT):
    '''Reconcile a single unrooted gene tree. (Iterates over all rooted versions of this).'''

    initFamNum,speciesRtreeO,geneUtreeO,tipMapD,gtLocusMapD,D,T,L,O,R = argT
    
    # convert species tree to dp format
    speciesTreeD = speciesRtreeO.createDtlorD(True)

    # try all rootings and record the ones and their
    # solutions with the best scores
    minCost=float('inf')
    bestRootingsL=[]
    for geneRtreeO in geneUtreeO.iterAllRootedTreesIncludeBranchLen():

        cost,G = reconcileOneRootedGeneTree(geneRtreeO,speciesTreeD,tipMapD,gtLocusMapD,D,T,L,O,R)
        
        if cost<minCost: 
            #if the score is better than current best
            #update best score, clear record and add new record
            minCost=cost
            bestRootingsL=[]
            bestRootingsL.append((geneRtreeO,G,minCost))
        elif cost==minCost:
            bestRootingsL.append((geneRtreeO,G,minCost))

    #sample one G from the Gs for this specific unrooted tree
    optGeneRtreeO,optG,minCost=random.choice(bestRootingsL) 

    return initFamNum,optGeneRtreeO,optG,minCost

def reconcileOneRootedGeneTree(geneRtreeO,speciesTreeD,tipMapD,gtLocusMapD,D,T,L,O,R):
    '''Reconcile a single rooted gene tree.'''
    geneTreeD = geneRtreeO.createDtlorD(False) # put in dp format
    cost, G = new_DTLOR_DP.compute_dtlor_graph(speciesTreeD,geneTreeD,tipMapD,gtLocusMapD,D,T,L,O,R)
    return cost,G

def reconcilePermissiveOrigin(paramD,initialFamiliesO,speciesRtreeO,locusMapD,genesO):
    '''Identifies initialFamilies where there is a tendency to insert
repeatedly in the same syntenic location. Re-runs reconciliation on
these families, with a new set of parameters which are permissive
to origin events.
    '''
    numProcesses = paramD['numProcesses']
    D=int(paramD["DTLRcostPermissiveOrigin"])
    T=int(paramD["DTLRcostPermissiveOrigin"])
    L=int(paramD["DTLRcostPermissiveOrigin"])
    O=int(paramD["originCostPermissiveOrigin"])
    R=int(paramD["DTLRcostPermissiveOrigin"])
    
    ## find families to redo
    iFamsToReconcileS = getIfamsToReconcilePermissiveOrigin(paramD,speciesRtreeO.leaves(),initialFamiliesO)

    ## Reconcile with permissive origin cost

    # geneTreeO in these families is rooted (having been run once
    # through dtlor). Unroot it. If only 1 tip, leave it as rooted
    for iFamO in iFamsToReconcileS:
        if type(iFamO.geneTreeO).__name__ == 'Rtree' and iFamO.geneTreeO.leafCount()>1:
            geneUtreeO = iFamO.geneTreeO.unroot()
            iFamO.geneTreeO = geneUtreeO
        
    initialFamiliesO = reconcileGeneTrees(iFamsToReconcileS,speciesRtreeO,initialFamiliesO,locusMapD,genesO,numProcesses,D,T,L,O,R)
    # make dtlorCost attribute negative, as indicator these were done permissively
    for iFamO in iFamsToReconcileS:
        if iFamO.dtlorCost != None:
            iFamO.dtlorCost = -iFamO.dtlorCost
    
    return initialFamiliesO

def getIfamsToReconcilePermissiveOrigin(paramD,strainNamesT,initialFamiliesO):
    '''Identify initial families that we should reconcile with costs
permissive to origin events. We do this by blasting our sequences
against a file of model seqs. Anything with similarity, should be done
with permissive-origin dtlor costs.
    '''

    reconcilePermissiveOriginGeneListPath = paramD['reconcilePermissiveOriginGeneListPath']

    # get geneToIfam dict
    gene2IfamD = {}
    for iFamO in initialFamiliesO.iterFamilies():
        for geneNum in iFamO.iterGenes():
            gene2IfamD[geneNum] = iFamO.famNum
            
    # identify ifams which these genes belong to
    with open(reconcilePermissiveOriginGeneListPath,'r') as f:
        S = f.read()
        geneL = S.rstrip().split()
    
    iFamsToReconcileS = set()
    for geneStr in geneL:
        geneNum = int(geneStr.split('_')[0])
        iFamsToReconcileS.add(initialFamiliesO.getFamily(gene2IfamD[geneNum]))

    return iFamsToReconcileS
    
#### Origin families

def createOriginFamiliesO(speciesRtreeO,initialFamiliesO,paramD,genesO):
    '''Create and return an originFamilies object, based on the initial families and recocniliations.'''

    originFamiliesO = Families(speciesRtreeO)

    for initFamO in initialFamiliesO.iterFamilies():

        if initFamO.geneCount()==1:
            # single gene, have no recon
            originFamiliesO = addOriginFamilyFromInitialFamiliesO(initFamO,False,originFamiliesO,genesO)
        elif len(initFamO.geneTreeO.multifurcatingNodes()) > 0:
            originFamiliesO = addOriginFamilyFromInitialFamiliesO(initFamO,True,originFamiliesO,genesO)
        else:
            initFamO,originFamiliesO = addOriginFamilyFromReconciliation(initFamO,originFamiliesO,paramD,genesO)
            initialFamiliesO.familiesD[initFamO.famNum] = initFamO # must stick back in here. LocusFams unchanged.
            
    return initialFamiliesO,originFamiliesO
        
def addOriginFamilyFromInitialFamiliesO(initFamO,isMultiFurc,originFamiliesO,genesO):
    '''For cases where there is no reconciliation, base origin family on
corresponding initial family.
    '''
    initFamNum = initFamO.famNum
    famNum=originFamiliesO.getNumFamilies() # num for new family
    if isMultiFurc:
        originFamiliesO.initializeFamily(famNum,initFamO.mrca,"origin",sourceFam=initFamNum)
    else:
        # single gene family
        geneRtreeO,reconD = createSingleGeneFamilyGeneTreeRecon(initFamO,genesO)
        originFamiliesO.initializeFamily(famNum,initFamO.mrca,"origin",geneTreeO=geneRtreeO,dtlorMprD=reconD,sourceFam=initFamNum)

    for initLocFamO in initFamO.getLocusFamilies():
        locFamNum = originFamiliesO.getNumLocusFamilies() # num for new loc fam
        lfO = LocusFamily(famNum,locFamNum,initLocFamO.lfMrca,initLocFamO.locusNum)
        for geneNum in initLocFamO.iterGenes():
            lfO.addGene(geneNum,genesO)
        originFamiliesO.addLocusFamily(lfO)

    return originFamiliesO

def createSingleGeneFamilyGeneTreeRecon(initFamO,genesO):
    '''Construct a gene tree and reconD for a single gene family.'''
    
    geneL=list(initFamO.iterGenes())
    assert(len(geneL)==1) # verify single gene
    geneNum=geneL[0]
    strainName = genesO.numToStrainName(geneNum)
    locusNum = initFamO.locusFamiliesL[0].locusNum

    # tree
    nodeConnectD = {str(geneNum):(ROOT_PARENT_NAME,)}
    geneRtreeO = Rtree(nodeConnectD,str(geneNum))

    # recon
    reconD = {}
    reconD[(geneNum,'b')] = [('O', strainName, 'b', locusNum)]
    reconD[(geneNum,'n')] = [('M', strainName, 'n', locusNum)]

    return geneRtreeO,reconD
    
def addOriginFamilyFromReconciliation(initFamO,originFamiliesO,paramD,genesO):
    '''Given a rooted gene tree and a reconciliation, add origin
families. One origin family for each origin event.
    '''

    # in the first pass at originFamily formation, dtlorMprD will be
    # none (and we chould arbirarily choose one). In the second pass,
    # it will have the MPR we want to use in it.

    if initFamO.dtlorMprD == None:
        # if initFamO has dtlorMprD equals None, then arbitrarily choose
        # a median MPR from the graph object.
        mprOrigFormatD,mprNodeFormatD = initFamO.getMprReconDFromGraph(originFamiliesO.speciesRtreeO.preorder(),paramD,True,False)
        initFamO.addMprD(mprOrigFormatD) # keep the mpr used for future reference
    else:
        # If it has a value in dtlorMprD, then use that MPR.
        mprNodeFormatD = initFamO.getMprReconDFromMpr(originFamiliesO.speciesRtreeO.preorder(),paramD)
        
    # add origin families according to this MPR
    productFamL = [] # to keep ofams that come from this ifam
    for speciesMrca,geneRtreeO,splitReconD in iterSplitOfamData(mprNodeFormatD,initFamO):
            
        famNum = originFamiliesO.getNumFamilies()
        locusFamNum = originFamiliesO.getNumLocusFamilies()
        originFamiliesO.initializeFamily(famNum,speciesMrca,"origin",geneTreeO=geneRtreeO,dtlorMprD=splitReconD,sourceFam=initFamO.famNum)
        
        for lfO in iterLocusFamiliesInOrigin(splitReconD,geneRtreeO,famNum,locusFamNum,genesO):
            originFamiliesO.addLocusFamily(lfO)
        
        productFamL.append(famNum)
        
    # update ifam
    initFamO.productFamT = tuple(productFamL)
    
    return initFamO,originFamiliesO

def iterSplitOfamData(mprNodeFormatD,initFamO):
    '''Helper that splits mprNodeFormatD according to origin events
yielding various bits of data needed to make an origin family.'''
    
    # get origin defined trees
    branchOriginL = getBranchesWithSpecifiedEvents(mprNodeFormatD,"O")
    originTreeL = splitTreeByOrigin(initFamO.geneTreeO,branchOriginL)
    
    # extract parts of mprNodeFormatD corresponding to each tree
    for geneRtreeO,splitReconD in getGeneTreeReconPairs(originTreeL,mprNodeFormatD):

        keyNB = (geneRtreeO.rootNode,'n')

        if geneRtreeO.isLeaf(geneRtreeO.rootNode):
            speciesMrca = geneRtreeO.rootNode
        else:
            speciesMrca = splitReconD[keyNB][0][1] # remember, the value in mprNodeFormatD is a list of events

        yield speciesMrca,geneRtreeO,splitReconD

def getBranchesWithSpecifiedEvents(mprNodeFormatD,event):
    '''Given a mprNodeFormatD keyed by (geneTreeLoc,geneTreeNB) tuples, extract
the branches (geneTreeLoc's) where an O even occurred.'''
    outL = []
    for key,value in mprNodeFormatD.items():
        geneTreeLoc,geneTreeNB = key
        eventsInEntryL = []
        for eventT in value:
            eventsInEntryL.append(eventT[0])
        if event in eventsInEntryL:
            outL.append(geneTreeLoc)
    return outL

def splitTreeByOrigin(geneRtreeO,branchOriginL):
    '''Take a gene tree, and a list of branches where origin events
occur. Split the tree into multiple trees, where each origin event
defines the root of a new tree. Note that if there is no origin event
at the root (ie they come later) this may mean that some parts of the
gene tree are not included.

    '''
    originTreeL = []
    for branch in branchOriginL:
        originTreeL.append(geneRtreeO.subtree(branch))
    return originTreeL

def getGeneTreeReconPairs(originTreeL,mprNodeFormatD):
    '''Given a list of geneTrees (new origin families), extract the
corresponding parts from the mprNodeFormatD for each.'''

    geneTreeReconPairL = []
    for geneRtreeO in originTreeL:
        splitReconD = {}
        for geneTreeLoc in geneRtreeO.preorder():
            nbKey = (geneTreeLoc,'b')
            if nbKey in mprNodeFormatD:
                splitReconD[nbKey] = mprNodeFormatD[nbKey]
            nbKey = (geneTreeLoc,'n')
            if nbKey in mprNodeFormatD:
                splitReconD[nbKey] = mprNodeFormatD[nbKey]
        geneTreeReconPairL.append((geneRtreeO,splitReconD))
    return geneTreeReconPairL
        
def iterLocusFamiliesInOrigin(splitReconD,geneRtreeO,famNum,locusFamNum,genesO):
    '''Iterator yielding locus families. Given a geneRtreeO for an origin
family, and the corresponding splitReconD, divide up further into
locusFamilies. Each locus family is defined by the most recent R event
(or O) in its history. famNum is the family number all these locus families will share. locusFamNum is the place to start the locus family numbers.
    '''
    branchRL = getBranchesWithSpecifiedEvents(splitReconD,"R")
    freeGeneL,locFamTL = splitTreeIntoLocusFamilies(geneRtreeO,geneRtreeO.rootNode,branchRL)
    if len(freeGeneL)>0:
        locFamTL.append((geneRtreeO.rootNode,tuple(freeGeneL))) # add last, O event can begin LF

    for locFamT in locFamTL:
        geneRtreeORBranch,locFamGenesT = locFamT
        lfReconRootKey = (geneRtreeORBranch,'b')
        # pull out the R or O event to get mrca
        _,lfSpeciesMrca,_,locusAtBottom = [val for val in splitReconD[lfReconRootKey] if val[0] in 'OR'][0]
        
        lfO=LocusFamily(famNum,locusFamNum,lfSpeciesMrca,int(locusAtBottom),lfReconRootKey)
        lfO.addGenes(locFamGenesT,genesO)
        locusFamNum += 1
        yield lfO


def splitTreeIntoLocusFamilies(geneRtreeO,node,branchRL):
    '''Split tips in geneRtreeO into locus families. Each gene is put into
a locus family according to the most recent R event in its lineage. If
no R events, then defined by its O event (each tree has an O at
root). We do this recursively. Function returns a tuple
(outFreeGeneL,outLocFamTL). outFreeGeneL is simply a list of genes
that have not yet encountered an R event. outLocFamTL is a list of
tuples, where each tuple is a locus family. It contains two
elements. The left is the branch on which the defining R event
occurs. The right is a list of genes defined by this. If we are at a
node where an R event occurred (is present in branchRL) then we take
the genes in outFreeGeneL, make into a tuple, and add to outLocFamTL.

    '''
    # in gene tree tip number is string. when making a list of genes,
    # we should make it into an int
    if geneRtreeO.isLeaf(node) and node in branchRL:
        return [],[(node,(int(node),))]
    elif geneRtreeO.isLeaf(node) and node not in branchRL:
        return [int(node)],[] 
    else:
        childFreeGeneL = []
        childLocFamTL = []
        for childNode in geneRtreeO.children(node):
            tempFreeGeneL,tempLocFamTL = splitTreeIntoLocusFamilies(geneRtreeO,childNode,branchRL)
            childFreeGeneL.extend(tempFreeGeneL)
            childLocFamTL.extend(tempLocFamTL)
        if node in branchRL:
            outLocFamT = (node,tuple(childFreeGeneL))
            return [],[outLocFamT]+childLocFamTL
        else:
            return childFreeGeneL,childLocFamTL

#### Refine families

def refineFamilies(paramD,islandByNodeD,initialFamiliesO,originFamiliesO,geneOrderD,genesO,outputSummaryF,strainNamesT):
    '''Refine origin families by considering alternate reconciliations. We
determine what the best MPR to use is, and update initialFamiliesO
with this. We then re-run origin family formation from scratch and
create a new origin families object. Returns both initialFamiliesO and
originFamiliesO.
    '''
    upperNumMprThreshold = paramD ['upperNumMprThreshold']
    islandLenThresholdRefineFamilies = paramD['islandLenThresholdRefineFamilies']
    geneProximityRangeRefineFamilies = paramD['geneProximityRangeRefineFamilies']
    proximityThreshold = paramD['proximityThresholdMerge'] # larger than initial
    rscThreshold = paramD['rscThresholdMerge']
    geneProximityRange = paramD['geneProximityRange'] # standard one for island merging
    initFamilyFN = paramD['initFamilyFN']
    originFamilyFN =  paramD['originFamilyFN']

    speciesRtreeO = originFamiliesO.speciesRtreeO
    
    ## get candidates

    # get all ofams from locus islands of size
    # refineFamIslandLenThreshold or smaller. Then get the ifams which
    # these come from. The ones among these with multiple MPRS we will
    # call refineCandidateIfams
    numIfams1_2_famIslands=0
    refineCandidateIfamS = set()
    for islandsAtNodeL in islandByNodeD.values():
        for locusIslandO in islandsAtNodeL:
            if len(locusIslandO) <= islandLenThresholdRefineFamilies:
                for origLfO in locusIslandO.iterLocusFamilies(originFamiliesO):
                    origFamO = originFamiliesO.getFamily(origLfO.famNum)
                    # get sourceIfam
                    sourceIfamO = initialFamiliesO.getFamily(origFamO.sourceFam)
                    # we can only consider those with reconciliations
                    # (single gene families lack them, and for now so
                    # do families with multifurcating nodes)
                    if sourceIfamO.dtlorGraphD != None:
                        numIfams1_2_famIslands+=1
                        if sourceIfamO.countMPRs() > 1:
                            refineCandidateIfamS.add(sourceIfamO)

    # create geneToOfamD
    geneToOfamD = {}
    for origFamO in originFamiliesO.iterFamilies():
        for geneNum in origFamO.iterGenes():
            geneToOfamD[geneNum] = origFamO.famNum
                            
    ## Refine
    geneProximityD = genomes.createGeneProximityD(geneOrderD,geneProximityRange)

    # make argument list
    argumentL = []
    for candIfamO in refineCandidateIfamS:
        argumentL.append((candIfamO,geneOrderD,geneProximityRangeRefineFamilies,geneToOfamD,originFamiliesO,upperNumMprThreshold,speciesRtreeO,paramD,genesO,geneProximityD,proximityThreshold,rscThreshold))

    # run on multiple processors
    with Pool(processes=paramD['numProcesses']) as p:
        for ifamNum,bestMprOrigFormatD in p.imap_unordered(getBestMprOrigFormatD, argumentL):
            # put output bestMprs back in ifam objects
            ifam = initialFamiliesO.getFamily(ifamNum)
            ifam.addMprD(bestMprOrigFormatD)

    ## make origin families again
    # the mprs in some families in initialFamiliesO have been changed, so we'll get a different result here
    initialFamiliesO,originFamiliesO = createOriginFamiliesO(speciesRtreeO,initialFamiliesO,paramD,genesO)

    ## update summary file
    print("Family refinement:",file=outputSummaryF)
    summaryL = []
    summaryL.append(["Number of origin families in small islands",str(numIfams1_2_famIslands)])
    summaryL.append(["Number of these with multiple MPRs",str(len(refineCandidateIfamS))])
    printTable(summaryL,indent=2,fileF=outputSummaryF)
    
    print("Origin families after refinement:",file=outputSummaryF)
    writeFamilyFormationSummary(originFamiliesO,outputSummaryF)

    # write familes to file
    writeFamilies(initialFamiliesO,initFamilyFN,genesO,strainNamesT,paramD) # update with different mprs
    writeFamilies(originFamiliesO,originFamilyFN,genesO,strainNamesT,paramD)
    #writeFamilies(initialFamiliesO,"ifamRefine.out",genesO,strainNamesT,paramD) # update with different mprs
    #writeFamilies(originFamiliesO,"ofamRefine.out",genesO,strainNamesT,paramD)

    return initialFamiliesO,originFamiliesO

def getBestMprOrigFormatD(argumentL):
    '''Given a candidate ifam, find the best MPR given nearbyOfamL.'''

    # Note: passing the largish objects in argumentL could become a problem with large datasets

    candIfamO,geneOrderD,geneProximityRangeRefineFamilies,geneToOfamD,originFamiliesO,upperNumMprThreshold,speciesRtreeO,paramD,genesO,geneProximityD,proximityThreshold,rscThreshold = argumentL

    nearbyOfamL = getNearbyOfamL(candIfamO,geneOrderD,geneProximityRangeRefineFamilies,geneToOfamD,originFamiliesO)

    bestMprOrigFormatD,bestOfamL = getBestOfamsFromCandIfam(candIfamO,upperNumMprThreshold,speciesRtreeO,paramD,max(originFamiliesO.familiesD.keys()),genesO,nearbyOfamL,geneProximityD,proximityThreshold,rscThreshold)

    return candIfamO.famNum,bestMprOrigFormatD
        
def getNearbyOfamL(candIfamO,geneOrderD,geneProximityRangeRefineFamilies,geneToOfamD,originFamiliesO):
    '''Given a list of target ofams whose placement we're reconsidering,
get a collection of other ofams that are nearby to these.'''

    # get target Ofams
    targetOfamS = set((originFamiliesO.getFamily(ofnum) for ofnum in candIfamO.productOfams()))
    
    # get genes in target ofams
    targetGenesS = set()
    targetOfamNumS = set()
    for tofamO in targetOfamS:
        targetGenesS.update(tofamO.getAllGenes())
        targetOfamNumS.add(tofamO.famNum)

    # get nearby ofams
    nearbyOfamNumS = set()
    for contigT in geneOrderD.values():
        for geneNumT in contigT:
            for i in range(len(geneNumT)):
                if geneNumT[i] in targetGenesS:
                    # found a target gene, slice out genes around it
                    for nearbyGene in geneNumT[i-geneProximityRangeRefineFamilies:i+geneProximityRangeRefineFamilies+1]:
                        ofamNum = geneToOfamD[nearbyGene]
                        nearbyOfamNumS.add(ofamNum)
    # remove target fams themselves
    nearbyWithoutTargetNumS = nearbyOfamNumS - targetOfamNumS
    nearbyWithoutTargetL = []
    for ofamNum in nearbyWithoutTargetNumS:
        nearbyWithoutTargetL.append(originFamiliesO.getFamily(ofamNum))

    return nearbyWithoutTargetL

def getBestOfamsFromCandIfam(candIfamO,upperNumMprThreshold,speciesRtreeO,paramD,maxOfamNum,genesO,nearbyOfamL,geneProximityD,proximityThreshold,rscThreshold):
    '''Given an inital families object with multiple MPRs, determine the
best MPR by running island formation with nearby ofams. In the case
that there are more MPRs than upperNumMprThreshold, we randomly sample
from the space of MPRs.
    '''
    # get place to start numbering locus fams (so no collisions)
    maxLocFamNum = 0
    for ofamO in nearbyOfamL:
        for lfO in ofamO.getLocusFamilies():
            if lfO.locusFamNum > maxLocFamNum:
                maxLocFamNum = lfO.locusFamNum
    
    bestNumIslands = float('inf')
    bestCandMprOfamL = []
    bestMprOrigFormatD = {}
    testS = set()
    for mprOrigFormatD,candMprOfamL in iterCandidateMprOfams(candIfamO,upperNumMprThreshold,speciesRtreeO,paramD,maxOfamNum,maxLocFamNum,genesO):

        testFamiliesO = createFamiliesOFromListOfFamilies(nearbyOfamL+candMprOfamL,speciesRtreeO)
        locIslByNodeD = islands.createLocIslByNodeD(testFamiliesO,speciesRtreeO)

        testAllLocIslandsL = []
        for node in locIslByNodeD:
            if locIslByNodeD[node] != []:
                subRtreeO = speciesRtreeO.subtree(node)
                argT = (locIslByNodeD[node],geneProximityD,proximityThreshold,rscThreshold,subRtreeO,testFamiliesO)
                testAllLocIslandsL.extend(islands.mergeLocIslandsAtNode(argT))

        if len(testAllLocIslandsL) < bestNumIslands:
            bestNumIslands = len(testAllLocIslandsL)
            bestCandMprOfamL = candMprOfamL
            bestMprOrigFormatD = mprOrigFormatD

    return bestMprOrigFormatD,bestCandMprOfamL

def iterCandidateMprOfams(candIfamO,upperNumMprThreshold,speciesRtreeO,paramD,maxOfamNum,maxLocFamNum,genesO):
    '''Given an ifam object, iterate through MPRs, yielding the origin
families associated with each MPR. maxOfamNum specifies where to
start the numbering for famNum in these ofams. (so as not to clash
with what is present in other families we will merge with).
    '''

    ofamNum = maxOfamNum
    locusFamNum = maxLocFamNum

    if candIfamO.countMPRs() < upperNumMprThreshold:
        # iterate though all. Do not restrict to median mprs
        for mprOrigFormatD,mprNodeFormatD in candIfamO.iterMprReconDFromGraph(speciesRtreeO.preorder(),paramD,False):
            candMprOfamL = getCandMprOfamL(mprNodeFormatD,candIfamO,ofamNum,locusFamNum,genesO)
            yield mprOrigFormatD,candMprOfamL
    else:
        for i in range(upperNumMprThreshold):
            # take upperNumMprThreshold samples
            # randomly choose, do not restrict to median mprs
            mprOrigFormatD,mprNodeFormatD = candIfamO.getMprReconDFromGraph(speciesRtreeO.preorder(),paramD,False,True)
            candMprOfamL = getCandMprOfamL(mprNodeFormatD,candIfamO,ofamNum,locusFamNum,genesO)
            yield mprOrigFormatD,candMprOfamL

def getCandMprOfamL(mprNodeFormatD,candIfamO,ofamNum,locusFamNum,genesO):
    '''Get origin families associated with a particular candidate MPR.'''
    candMprOfamL = []
    for speciesMrca,geneRtreeO,splitReconD in iterSplitOfamData(mprNodeFormatD,candIfamO):
        ofamO = originFamily(ofamNum,speciesMrca,geneTreeO=geneRtreeO,dtlorMprD=splitReconD,sourceFam=candIfamO.famNum)
        # add locus families
        for lfO in iterLocusFamiliesInOrigin(splitReconD,geneRtreeO,ofamNum,locusFamNum,genesO):
            ofamO.addLocusFamily(lfO)
            locusFamNum += 1

        candMprOfamL.append(ofamO)

        ofamNum += 1
    return candMprOfamL
            
def createFamiliesOFromListOfFamilies(ofamL,speciesRtreeO):
    '''Create a familiesO object with the origin families in the input list.'''
    testFamiliesO = Families(speciesRtreeO)
    for ofamO in ofamL:
        # add family
        testFamiliesO.initializeFamily(ofamO.famNum,ofamO.mrca,"origin",geneTreeO=ofamO.geneTreeO,dtlorMprD=ofamO.dtlorMprD,sourceFam=ofamO.sourceFam)
        # add locus families
        for lfO in ofamO.getLocusFamilies():
            testFamiliesO.addLocusFamily(lfO)
    return testFamiliesO
    
#### Input/output

def writeFamilies(familiesO,familyFN,genesO,strainNamesT,paramD):
    '''Write all gene families to familyFN, one family per line.'''

    geneInfoFN = paramD['geneInfoFN']

    # get the num to name dict, only for strains we're looking at.
    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT))
    
    with open(familyFN,'w') as f:
        for fam in familiesO.iterFamilies():
            f.write(fam.fileStr(genesO)+'\n')

def readFamilies(familyFN,speciesRtreeO,genesO,famType):
    '''Read the family file named familyFN, creating a Families object.
    '''
    
    familiesO = Families(speciesRtreeO)
    f=open(familyFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.split('\t')
        famNum=int(L[0])
        mrca = L[1]

        # tree
        if L[2] == "None":
            geneRtreeO = None
        else:
            geneRtreeO = Rtree() # currently always rooted
            geneRtreeO.fromString(L[2])

        # dtlor cost
        if L[3] == "None":
            dtlorCost = None
        else:
            dtlorCost = int(L[3])

        # dtlorGraphD
        if L[4] == "None":
            dtlorGraphD = None
        else:
            dtlorGraphD = eval(L[4])

        # dtlorMprD
        if L[5] == "None":
            dtlorMprD = None
        else:
            dtlorMprD = eval(L[5])
        
        # sourceFam
        if L[6] == "None":
            sourceFam = None
        else:
            sourceFam = int(L[6])

        # productFamT
        if L[7] == "None":
            productFamT = None
        else:
            productFamT = eval(L[7])

        # make it
        familiesO.initializeFamily(famNum,mrca,famType,geneTreeO=geneRtreeO,dtlorCost=dtlorCost,dtlorGraphD=dtlorGraphD,dtlorMprD=dtlorMprD,sourceFam=sourceFam,productFamT=productFamT)

        # add locus families
        lfL = L[8:]
        for lfStr in lfL:
            lfSplitL = lfStr.rstrip().split(',')
            locusFamNum=int(lfSplitL[0])
            lfMrca = lfSplitL[1]
            locusNum = eval(lfSplitL[2])
            
            if lfSplitL[3] == 'None':
                reconRootKey = None
            else:
                geneTreeLoc,geneTreeNB = lfSplitL[3].split('_')
                reconRootKey = (geneTreeLoc,geneTreeNB)
            
            geneL=[]
            for geneName in lfSplitL[4:]:
                geneNum = int(geneName.split('_')[0])
                geneL.append(geneNum)
            lfO = LocusFamily(famNum,locusFamNum,lfMrca,locusNum,reconRootKey)
            lfO.addGenes(geneL,genesO)
            familiesO.addLocusFamily(lfO)

    f.close()
    return familiesO
