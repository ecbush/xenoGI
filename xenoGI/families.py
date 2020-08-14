import sys,numpy,os,random,glob,copy,shutil
sys.setrecursionlimit(100000)
from scipy.signal import find_peaks
from Bio import Phylo
from multiprocessing import Pool
from . import trees,scores,DTLOR_DP,new_DTLOR_DP
from .Family import *
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

    # threshold for maximum size of gene tree, if above, split
    maxIfamSize = paramD['maxIfamSize']

    # define some variables
    initFamilyFN = paramD['initFamilyFN']
    originFamilyFN =  paramD['originFamilyFN']
    geneInfoFN = paramD['geneInfoFN']

    # checks
    homologyCheck(genesO,aabrhHardCoreL,scoresO,outputSummaryF,paramD)

    # create blast families, output is directory of gene trees
    createBlastFamilies(paramD,scoresO,genesO,outputSummaryF)
    print("Finished making gene trees",file=sys.stderr)

    initialFamiliesO,locusMapD = createInitialFamiliesO(paramD,genesO,aabrhHardCoreL,scoresO,speciesRtreeO,outputSummaryF)
    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT)) # needed for writeFamilies
    writeFamilies(initialFamiliesO,initFamilyFN,genesO,strainNamesT,paramD)
    print("Initial families:",file=outputSummaryF)
    writeFamilyFormationSummary(initialFamiliesO,outputSummaryF)
    
    # reconcile
    initialFamiliesO = reconcileAllGeneTrees(speciesRtreeO,initialFamiliesO,locusMapD,genesO,paramD)
    writeFamilies(initialFamiliesO,initFamilyFN,genesO,strainNamesT,paramD)
    print("Initial families updated with reconciliations",file=sys.stderr)    

    # create origin families
    originFamiliesO = createOriginFamiliesO(speciesRtreeO,initialFamiliesO,paramD,genesO)

    # write origin familes to file
    writeFamilies(originFamiliesO,originFamilyFN,genesO,strainNamesT,paramD)
    print("Origin families:",file=outputSummaryF)
    writeFamilyFormationSummary(originFamiliesO,outputSummaryF)

    return originFamiliesO

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

#### Create initial families object

def createBlastFamilies(paramD,scoresO,genesO,outputSummaryF):
    '''Given a scoresO object, create gene families based on blast
    connectivity, then make gene trees with these.

    '''
    geneFamilyTreesDir = paramD['geneFamilyTreesDir']
    blastFamGeneTreeFileStem = paramD['blastFamGeneTreeFileStem']
    blastFamilyFN = paramD['blastFamilyFN']
    
    ## get blast families as list of sets
    blastFamilySetL = createBlastFamilySetL(scoresO, genesO)
    print("  Number of blast families:",len(blastFamilySetL),file=outputSummaryF)
    blastFamilySetL.sort(reverse=True,key=len) # sort by number of genes in descending order
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

    trees.makeGeneTrees(paramD,False,genesO,geneFamilyTreesDir,blastFamGeneTreeFileStem,blastFamilySetL)

    # remove alignments
    for fn in glob.glob(os.path.join(geneFamilyTreesDir,"align*.afa")):
        os.remove(fn)
    
def createBlastFamilySetL(scoresO, genesO):
    '''
    Input
    ------------------------------------------------
    scoresO:        a Score object that scores the 
                    BLAST results from all the genes

    Output
    ------------------------------------------------
    initFamilySetL:      a list of sets where each set stores
                    all the genes that are connected as
                    indicated by significant BLAST score 
                    (that appear in scoresO)
    '''
   
    def stronglyConnected(temp, gene, visited): 
        visited.add(gene)
        temp.add(gene) 
  
        # Repeat for all vertices adjacent 
        # to this gene
        neighbors=scoresO.getConnectionsGene(gene)
        degreeD[gene]=len(neighbors)
        if neighbors:
            for i in neighbors: 
                if i not in visited:     
                    # Update the list 
                    temp = stronglyConnected(temp, i, visited) 
        
        return temp

    # main for createBlastFamilySetL
    degreeD = {}
    allGenesL=list(genesO.iterGenes())
    scoresO.createNodeConnectD() 
    connectedGenes=list(scoresO.nodeConnectD.keys())
    blastFamilySetL=[]
    visited=set()
    for gene in allGenesL: 
        if gene in scoresO.nodeConnectD:
            if gene not in visited: 
                temp =set()
                newFam=stronglyConnected(temp, gene, visited)
                
                blastFamilySetL.append(newFam) 
        else:
            fam=set()
            fam.add(gene)
            blastFamilySetL.append(fam)
  
    return blastFamilySetL

def createInitialFamiliesO(paramD,genesO,aabrhHardCoreL,scoresO,speciesRtreeO,outputSummaryF):
    ''''''
    blastFamGeneTreeFileStem = paramD['blastFamGeneTreeFileStem']
    aabrhHardCoreGeneTreeFileStem = paramD['aabrhHardCoreGeneTreeFileStem']
    geneFamilyTreesDir = paramD['geneFamilyTreesDir']
    maxIfamSize = paramD['maxIfamSize']
    forceSplitUtreeBalanceMultiplier = paramD['forceSplitUtreeBalanceMultiplier']
    
    # load gene trees
    blastFamGeneTreeL = loadGeneTrees(paramD,blastFamGeneTreeFileStem)

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
    trees.makeGeneTrees(paramD,False,genesO,geneFamilyTreesDir,aabrhHardCoreGeneTreeFileStem,newAabrhHardCoreL)

    aabrhHardCoreGeneTreeL = loadGeneTrees(paramD,aabrhHardCoreGeneTreeFileStem)

    # remove alignments
    for fn in glob.glob(os.path.join(geneFamilyTreesDir,"align*.afa")):
        os.remove(fn)

    # obtain threshold for utree splitting
    splitThresh = calculateTreeSplitThreshold(paramD,aabrhHardCoreGeneTreeL)

    # split on branches larger than splitThresh
    singleGeneTreeL = [] # need this since splitting funcs don't like single tip trees
    geneTreeL = []
    for blastFamNum,blastFamUtreeO in blastFamGeneTreeL:
        if blastFamUtreeO.leafCount() == 1:
            singleGeneTreeL.append((blastFamUtreeO.leafCount(),blastFamNum,blastFamUtreeO))
        else:
            for aUtO in splitUtreeThreshold([blastFamUtreeO],splitThresh)           :
                if aUtO.leafCount() == 1:
                    singleGeneTreeL.append((aUtO.leafCount(),blastFamNum,aUtO))
                else:
                    # failsafe to keep below maxIfamSize. If none too
                    # big, will pass back same list
                    for bUtO in splitUtreeFailsafe([aUtO],maxIfamSize,forceSplitUtreeBalanceMultiplier):
                        geneTreeL.append((bUtO.leafCount(),blastFamNum,bUtO))

    geneTreeL.extend(singleGeneTreeL)
    del singleGeneTreeL
    geneTreeL.sort(reverse=True,key=lambda x: x[0]) # big families first
    
    # make initial families object
    
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

def loadGeneTrees(paramD,geneTreeFileStem):
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

def splitUtreeThreshold(utreeL,splitThresh):
    '''Given an input list containing Utree objects, recursively split
these until there are no branches longer than splitThresh. Trees in
utreeL should not be single tip trees.
    '''

    def getIndFirstAboveThresh(utreeL):
        i=-1 # for empty list case
        for i in range(len(utreeL)):
            maxBranchPair,maxBrLen = utreeL[i].maxBranchLen()
            if maxBrLen > splitThresh:
                return i
        return i+1 # none above thresh
        
    # main splitUtree
    aboveThreshInd = getIndFirstAboveThresh(utreeL)

    if aboveThreshInd == len(utreeL): # done
        return utreeL
    else:
        completeL = utreeL[:aboveThreshInd] # these are done
        
        # split the branch we found that was too large
        maxBranchPair,maxBrLen = utreeL[aboveThreshInd].maxBranchLen()

        remainderL = []
        aUtreeO,bUtreeO = utreeL[aboveThreshInd].split(maxBranchPair)

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
            
        # put the results back in the todo list and recurse
        remainderL.extend(utreeL[aboveThreshInd+1:])
        remainderSplitL = splitUtreeThreshold(remainderL,splitThresh)
        completeL.extend(remainderSplitL)
        return completeL

def splitUtreeFailsafe(utreeL,maxIfamSize,forceSplitUtreeBalanceMultiplier):
    '''Function for cutting tree size down in order to limit dtlor
calculation time.  Given an input list containing Utree objects,
recursively split these until all are below maxIfamSize. Trees in
utreeL should not be single tip trees.
    '''

    def getIndFirstTooBig(utreeL):
        i=-1 # for empty list case
        for i in range(len(utreeL)):
            if utreeL[i].leafCount() > maxIfamSize:
                return i
        return i+1 # none too big
        
    # main splitUtreeFailsafe
    tooBigInd = getIndFirstTooBig(utreeL)

    if tooBigInd == len(utreeL): # done
        return utreeL
    else:
        completeL = utreeL[:tooBigInd] # these are done
        
        # split the tree that was too large
        aUtreeO,bUtreeO = forceSplitUtree(utreeL[tooBigInd],forceSplitUtreeBalanceMultiplier)
        
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
        remainderSplitL = splitUtreeFailsafe(remainderL,maxIfamSize,forceSplitUtreeBalanceMultiplier)
        completeL.extend(remainderSplitL)
        return completeL

def forceSplitUtree(geneUtreeO,forceSplitUtreeBalanceMultiplier):
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
    branchLenL = list(geneUtreeO.branchLenD.items())
    branchLenL.sort(reverse=True,key=lambda x: x[1])
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

    # split
    splitBranch = comboL[0][0]
    #print(splitBranch,geneUtreeO.branchLenD[splitBranch])
    aO,bO = geneUtreeO.split(splitBranch)
    return aO,bO

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
    initialFamiliesO.initializeFamily(famNumCounter,mrca,"initial",geneUtreeO,None,sourceFam)
    for locusFamilyL in locusFamLL: #locusFamilyL contains all the genes in that lf
        speciesL=[]
        for gene in locusFamilyL:
            locusMapD[gene]=locusFamNumCounter
            speciesL.append(genesO.numToStrainName(gene))
        lfMrca=speciesRtreeO.findMrca(speciesL)
        lf=LocusFamily(famNumCounter,locusFamNumCounter,lfMrca)
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
    def stronglyConnected(tempFamL, gene, visitedS,geneNeighborsD):
        '''Find components connected to gene and put in tempFamL.'''
        visitedS.add(gene)
        tempFamL.append(gene)
        neighborsL=geneNeighborsD[gene]
        if len(neighborsL)>0:
            for neighbGene in neighborsL: 
                if neighbGene not in visitedS:     
                    # Update the list 
                    tempFamL = stronglyConnected(tempFamL, neighbGene, visitedS, geneNeighborsD) 
        return tempFamL 

    locusFamLL=[]
    visitedS=set()
    for gene in genesL: 
        if gene not in visitedS: 
            tempFamL =[]
            newFamL=stronglyConnected(tempFamL, gene, visitedS, geneNeighborsD)
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
    for fam in familiesO.iterFamilies():
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

    summaryL.append(["Number of families with one LocusFamily",str(singleLfFams)])
    summaryL.append(["Number of families with multiple LocusFamilies",str(multipleLfFams)])
    summaryL.append(["Number of families with gene(s) in only one strain",str(singleStrainFams)])
    summaryL.append(["Number of families with genes in multiple strains",str(multipleStrainFams)])
    
    printTable(summaryL,indent=2,fileF=outputSummaryF)

    return

#### Reconciliation

def reconcileAllGeneTrees(speciesRtreeO,initialFamiliesO,locusMapD,genesO,paramD):
    '''Reconcile gene family trees to the species tree using the DTLOR algorithm.'''

    D=int(paramD["duplicationCost"])
    T=int(paramD["transferCost"])
    L=int(paramD["lossCost"])
    O=int(paramD["originCost"])
    R=int(paramD["rearrangeCost"])

    argumentL = []
    for ifam in initialFamiliesO.iterFamilies():
        initFamNum = ifam.famNum
        geneUtreeO = ifam.geneTreeO

        assert(ifam.geneCount() == geneUtreeO.leafCount()) # temp

        # skip single gene and multifurcating families
        if ifam.geneCount() == 1 or len(geneUtreeO.multifurcatingNodes()) > 0:
            continue
        
        # in loop
        tipMapD=getTipMapping(geneUtreeO,genesO)

        # Make new gtLocusMapD with only those genes in geneUtreeO
        gtLocusMapD = reduceLocusMap(geneUtreeO,locusMapD)
        
        # add to argumentL
        argT = (initFamNum,speciesRtreeO,geneUtreeO,tipMapD,gtLocusMapD,D,T,L,O,R)
        argumentL.append(argT)

    # run on multiple processors
    with Pool(processes=paramD['numProcesses']) as p:
        for initFamNum,optGeneRtreeO,optG,minCost in p.imap_unordered(reconcile, argumentL):
            
            # store
            ifam = initialFamiliesO.getFamily(initFamNum)
            ifam.addGeneTree(optGeneRtreeO)
            ifam.addReconciliation(optG)
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
        
def reconcile(argT):
    '''Reconcile a single gene tree.'''

    initFamNum,speciesRtreeO,geneUtreeO,tipMapD,gtLocusMapD,D,T,L,O,R = argT
    
    # convert species tree to dp format
    speciesTreeD = speciesRtreeO.createDtlorD(True)

    # try all rootings and record the ones and their
    # solutions with the best scores
    minCost=float('inf')
    bestRootingsL=[]
    for geneRtreeO in geneUtreeO.iterAllRootedTrees():

        geneTreeD = geneRtreeO.createDtlorD(False) # put in dp format
        cost, G = new_DTLOR_DP.compute_dtlor_graph(speciesTreeD,geneTreeD,tipMapD,gtLocusMapD,D,T,L,O,R)

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

    #MPR = new_DTLOR_DP.find_MPR(G)
    #print("G,MPR",sys.getsizeof(G),sys.getsizeof(MPR))

    return initFamNum,optGeneRtreeO,optG,minCost
    
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
            # add in families from reconciliation
            sourceFam = initFamO.famNum
             # get arbitrarily chosen median mpr
            reconD = initFamO.getMedianMprReconD(speciesRtreeO.preorder(),paramD,False)
            originFamiliesO = addOriginFamilyFromReconciliation(initFamO.geneTreeO,reconD,originFamiliesO,sourceFam,genesO)

    return originFamiliesO
        
def addOriginFamilyFromInitialFamiliesO(initFamO,isMultiFurc,originFamiliesO,genesO):
    '''For cases where there is no reconciliation, base origin family on
corresponding initial family.
    '''
    initFamNum = initFamO.famNum
    famNum=originFamiliesO.getNumFamilies() # num for new family
    if isMultiFurc:
        originFamiliesO.initializeFamily(famNum,initFamO.mrca,"origin")
    else:
        # single gene family
        geneRtreeO,reconD = createSingleGeneFamilyGeneTreeRecon(initFamO,genesO)
        originFamiliesO.initializeFamily(famNum,initFamO.mrca,"origin",geneRtreeO,reconD)

    originFamiliesO.getFamily(famNum).sourceFam = initFamNum

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
    
def addOriginFamilyFromReconciliation(geneRtreeO,reconD,originFamiliesO,sourceFam,genesO):
    '''Given a rooted gene tree and a reconciliation, add origin
families. One origin family for each origin event.

    '''
    # get origin defined trees
    branchOriginL = getBranchesWithSpecifiedEvents(reconD,"O")
    originTreeL = splitTreeByOrigin(geneRtreeO,branchOriginL)
    
    # extract parts of reconD corresponding to each tree
    for geneRtreeO,splitReconD in getGeneTreeReconPairs(originTreeL,reconD):
        famNum = originFamiliesO.getNumFamilies()
        keyNB = (geneRtreeO.rootNode,'n')

        if geneRtreeO.isLeaf(geneRtreeO.rootNode):
            speciesMrca = geneRtreeO.rootNode
        else:
            speciesMrca = splitReconD[keyNB][0][1] # remember, the value in reconD is a list of events
            
        originFamiliesO.initializeFamily(famNum,speciesMrca,"origin",geneRtreeO,splitReconD)
        originFamiliesO.getFamily(famNum).sourceFam = sourceFam
        originFamiliesO = getLocusFamiliesInOrigin(splitReconD,geneRtreeO,originFamiliesO,famNum,genesO)

    return originFamiliesO

def getBranchesWithSpecifiedEvents(reconD,event):
    '''Given a reconD keyed by (geneTreeLoc,geneTreeNB) tuples, extract
the branches (geneTreeLoc's) where an O even occurred.'''
    outL = []
    for key,value in reconD.items():
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

def getGeneTreeReconPairs(originTreeL,reconD):
    '''Given a list of geneTrees (new origin families), extract the
corresponding parts from the reconD for each.'''

    geneTreeReconPairL = []
    for geneRtreeO in originTreeL:
        splitReconD = {}
        for geneTreeLoc in geneRtreeO.preorder():
            nbKey = (geneTreeLoc,'b')
            if nbKey in reconD:
                splitReconD[nbKey] = reconD[nbKey]
            nbKey = (geneTreeLoc,'n')
            if nbKey in reconD:
                splitReconD[nbKey] = reconD[nbKey]
        geneTreeReconPairL.append((geneRtreeO,splitReconD))
    return geneTreeReconPairL
        
def getLocusFamiliesInOrigin(splitReconD,geneRtreeO,originFamiliesO,famNum,genesO):
    '''Given a geneRtreeO for an origin family, and the corresponding
splitReconD, divide up further into locusFamilies. Each locus family
is defined by the most recent R event (or O) in its history.

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
        newLocusNum=originFamiliesO.getNumLocusFamilies()
        lf=LocusFamily(famNum,newLocusNum,lfSpeciesMrca,int(locusAtBottom),lfReconRootKey)
        lf.addGenes(locFamGenesT,genesO)
        originFamiliesO.addLocusFamily(lf)
    return originFamiliesO

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
        
#### Input/output

def writeFamilies(familiesO,familyFN,genesO,strainNamesT,paramD):
    '''Write all gene families to fileName, one family per line.'''

    geneInfoFN = paramD['geneInfoFN']

    # get the num to name dict, only for strains we're looking at.
    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT))
    
    f=open(familyFN,'w')
    for fam in familiesO.iterFamilies():
        f.write(fam.fileStr(genesO)+'\n')
    f.close()

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

        if L[2] == "None":
            geneRtreeO = None
        else:
            geneRtreeO = Rtree()
            geneRtreeO.fromString(L[2])
        
        recon = eval(L[3].replace('inf', "float('inf')")) # a hack!
        # eval didn't like the string inf.
        sourceFam = eval(L[4])
        familiesO.initializeFamily(famNum,mrca,famType,geneRtreeO,recon,sourceFam)
        
        lfL = L[5:]
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
