# Functions for a modified version of the PhiGs algorithm
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# (We've added the use of synteny)
import sys,numpy
from scipy.signal import find_peaks
from . import trees,scores
from .Family import *
from .analysis import printTable

#### Main function

def createFamiliesO(tree,strainNum2StrD,scoresO,geneNames,aabrhL,paramD,subtreeL,outputSummaryF):
    '''Given a graph of genes and their similarity scores find families
using a PhiGs-like algorithm, with synteny also considered.

Here's a hierarchy of (major) function calls below this:
createFamiliesO
  createAllFamiliesDescendingFromInternalNode
    createFamilyFromSeed
  createAllFamiliesAtTip
  createAllLocusFamiliesDescendingFromInternalNode
     createAllLocusFamiliesAtOneInternalNode
        createLocusFamilyFromSeed
     createAllLocusFamiliesAtTip
    '''
    # Do an initial check the placement of the rootFocalClade on the
    # tree. If there are no outgroups, write a warning to our output
    # file
    if strainNum2StrD[tree[0]] == paramD['rootFocalClade']:
        print("""Warning: the chosen rootFocalClade falls at the root of the input
tree and thus does not have any outgroups. This is not recommended
because it can lead to problems accurately recognizing core gene
families in the presence of gene deletion."""+"\n",file=outputSummaryF)
    
    # initialize scoresO.nodeConnectL and scoresO.ScoreSummaryD for
    # use below
    strainNumsL=sorted(leaf for leaf in trees.leafList(tree))
    scoresO.createNodeConnectL(geneNames)
    scoresO.createAabrhScoreSummaryD(strainNumsL,aabrhL,geneNames)

    # create an object of class Families to store this in.
    familiesO = Families(tree)
    famNumCounter = 0
    locusFamNumCounter = 0

    # other assorted things we'll need
    geneUsedL = [False for x in geneNames.nums]
    nodeGenesL = createNodeGenesL(tree,geneNames) # has genes divided by node
    tipFamilyRawThresholdD = getTipFamilyRawThresholdD(tree,scoresO,paramD)


    # get homology thresholds for each strain pair in family
    # formation. Below this threshold we won't add a gene to
    homologousPeakMissingL,absMinRawThresholdForHomologyD = getAbsMinRawThresholdForHomologyD(paramD,scoresO)

    if homologousPeakMissingL != []:
        print("""Warning: for one or more strain pairs, we failed to find a homologous
(right) peak in the raw score histogram. For these cases we will use
default values for score thresholds in family formation. The pair(s)
in question are:""",file=outputSummaryF)
        for pairT in homologousPeakMissingL:
            print("  ",strainNum2StrD[pairT[0]]+'-'+strainNum2StrD[pairT[1]],file=outputSummaryF)
        
        print("""A possible explanation for this failure is that one or more species
is too distantly related. If this is the case it will result in poor
family formation, and most families (e.g. 80% or more) will have genes
in only one strain."""+"\n",file=outputSummaryF)

    # Create synThresholdD
    synThresholdD = getSynThresholdD(scoresO,tree,paramD)
        
    for familyMrca,lchild,rchild in createNodeProcessOrderList(tree):
        # this is preorder, so we get internal nodes before tips
        if lchild != None:
            # not a tip
                        
            geneUsedL,locusFamNumCounter,famNumCounter,familiesO = createAllFamiliesDescendingFromInternalNode(subtreeL,familyMrca,nodeGenesL,scoresO,geneNames,absMinRawThresholdForHomologyD,synThresholdD,paramD,geneUsedL,familiesO,famNumCounter,locusFamNumCounter)

        else:
            geneUsedL,locusFamNumCounter,famNumCounter,familiesO = createAllFamiliesAtTip(nodeGenesL,familyMrca,geneUsedL,tipFamilyRawThresholdD,scoresO,geneNames,absMinRawThresholdForHomologyD,synThresholdD,paramD,familiesO,famNumCounter,locusFamNumCounter)

       
    # Write family formation summary file

    summaryL=[]
    summaryL.append(["Total number of Families",str(len(familiesO.familiesD))])
    summaryL.append(["Total number of LocusFamilies",str(len(familiesO.locusFamiliesD))])

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
    
    printTable(summaryL,indent=0,fileF=outputSummaryF)
    
    writeFamilies(familiesO,geneNames,strainNum2StrD,paramD['familyFN'])

    return familiesO

## Support functions

def createNodeGenesL(tree,geneNames):
    '''Create a data structure to organize genes by strain. Returns a list
where the indices correspond to strain number and the elements are
sets.
    '''
    nodeGenesL=[set() for i in range(trees.nodeCount(tree))]
    for gene in geneNames.nums:
        strain = geneNames.numToStrainNum(gene)
        nodeGenesL[strain].add(gene)
    return nodeGenesL

def createNodeProcessOrderList(tree):
    '''Given a tree, output a list specifying nodes in pre-order (ancestors
before their descendants). For each node we give a tuple (node #, left
node #, right node #).
    '''
    if tree[1] == ():
        return [(tree[0], None, None)]
    else:
        l = createNodeProcessOrderList(tree[1])
        r = createNodeProcessOrderList(tree[2])
        return [(tree[0], tree[1][0], tree[2][0])] + l + r

def getTipFamilyRawThresholdD(tree,scoresO,paramD):
    '''Return a dictionary containing a raw score threshold for each
tip. This threshold is for use in forming families on the tip,
defining the minimum distance within which we will combine two genes
into a family.'''

    tipFamilyRawThresholdD = {}
    for leaf in trees.leafList(tree):
        # get average score at core genes for neighbors

        # put in call to get average...
        threshold,std = getNearestNeighborAverageScore(leaf,tree,scoresO)
        
        # multiply in an adjustment parameter (since average core gene
        # scores of neighbors would be way too high)
        threshold *= paramD['singleStrainFamilyThresholdAdjust']

        tipFamilyRawThresholdD[leaf] = threshold

    return tipFamilyRawThresholdD

def getNearestNeighborAverageScore(species,tree,scoresO):
    '''Get all the nearest neighbors of species, and collect the average
score for each against species at aabrh core genes. Return the average
of this. Assumes that scoreSummaryD has been initialized in
scoresO.
    '''
    neighbL = trees.getNearestNeighborL(species,tree)
    avScore = 0
    avStd = 0
    for neighb in neighbL:
        sc,std = scoresO.scoreSummaryD[(species,neighb)]
        avScore += sc
        avStd += std
    avScore /= len(neighbL)
    avStd /= len(neighbL)
    return avScore,avStd

## Histograms and thresholds

def getAbsMinRawThresholdForHomologyD(paramD,scoresO):
    '''For each pair of strains (including vs. self) determine a minimum
raw score below which we take scores to indicate non-homology. Return
in a dictionary keyed by strain pair. It also returns a boolean flag,
homologousPeakMissing. This will be true if any of the strain pairs
failed to yield a peak on the right side of the histogram
corresponding to homologous scores. This function works as
follows. For a given strain pair, it expects the histogram of scores
to have two peaks. A left non homologous peak, and a right homologous
peak. It attempts to find the right peak, and then finds the left
peak. It takes the right base of the left non-homologous peak to be
the threshold. Unless the left base of the right peak is actually
smaller, in which case it takes that.

    '''

    scoreHistNumBins = paramD['scoreHistNumBins']
    binWidth = 1.0/scoreHistNumBins # since scores range from 0-1

    homologousPeakMissingL = []
    homologyRawThresholdD = {}
    for strainPair in scoresO.getStrainPairs():

        scoreIterator = scoresO.iterateScoreByStrainPair(strainPair,'rawSc')
        binHeightL,indexToBinCenterL = scoreHist(scoreIterator,scoreHistNumBins)
        homologPeakLeftExtremePos=homologPeakChecker(binHeightL,indexToBinCenterL,binWidth,paramD)

        if homologPeakLeftExtremePos == float('inf'):
            homologousPeakMissingL.append(strainPair)
            # if we can't find the homologous peak, there's really no
            # point in examining the non-homologous peak. Just use
            # default threshold.
            threshold = paramD['defaultAbsMinRawThresholdForHomology']
        else:
            threshold = getMinRawThreshold(binHeightL,indexToBinCenterL,binWidth,homologPeakLeftExtremePos,paramD)
        homologyRawThresholdD[strainPair] = threshold
        homologyRawThresholdD[(strainPair[1],strainPair[0])] = threshold # flip order

    return homologousPeakMissingL,homologyRawThresholdD
    
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

def getMinRawThreshold(binHeightL,indexToBinCenterL,binWidth,homologPeakLeftExtremePos,paramD):
    '''Given a list of bin heights and another list giving the score
values at the middle of each bin, determine a threshold below which a
score should be taken to indicate non-homology. We do this by looking
for the non-homologous (left) peak in the score histogram. This
function assumes there is a right homologous peak present and takes
the left extreme of this peak as input.
    '''

    L = findPeaksOneCase(binHeightL,indexToBinCenterL,binWidth,paramD['nonHomologPeakWidth'],paramD['widthRelHeight'],paramD['nonHomologPeakProminence'],paramD['nonHomologLeftPeakLimit'],paramD['nonHomologRightPeakLimit'])

    if L == []:
        # no peak found, use default threshold
        threshold = paramD['defaultAbsMinRawThresholdForHomology']
    else:
        L.sort(reverse=True) # in the unlikely case there's more than one
        peakHeight,peakPos,leftExtremeOfPeakPos,rightExtremeOfPeakPos = L[0]

        # we now find the minimum of these two. We're after a threshold
        # where we're confident that things below it are definitely not
        # homologous.
        threshold = min(rightExtremeOfPeakPos,homologPeakLeftExtremePos)

    return threshold

def getSynThresholdD(scoresO,tree,paramD):
    '''Creates a dictionary to store synteny thresholds. This dictionary
itself contains two dictionaries one each for the minSynThresh
(minimum synteny score allowed for family formation) and the
synAdjustThreshold (the synteny score above which we adjust up the raw
score to make family formation more likely.) These dictionaries in
turn are keyed by strain pair.'''

    synThresholdD = {}
    synThresholdD['minSynThreshold'] = {}
    synThresholdD['synAdjustThreshold'] = {}
    
    for strainPair in scoresO.getStrainPairs():
        scoreIterator = scoresO.iterateScoreByStrainPair(strainPair,'synSc')
        binHeightL,indexToBinCenterL = scoreHist(scoreIterator,paramD['scoreHistNumBins'])

        minSynThreshold,synAdjustThreshold = getSynThresholds(binHeightL,indexToBinCenterL,paramD)

        synThresholdD['minSynThreshold'][strainPair] = minSynThreshold
        synThresholdD['minSynThreshold'][(strainPair[1],strainPair[0])] = minSynThreshold
        
        synThresholdD['synAdjustThreshold'][strainPair] = synAdjustThreshold
        synThresholdD['synAdjustThreshold'][(strainPair[1],strainPair[0])] = synAdjustThreshold

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
            
            minSynThreshold,synAdjustThreshold = getTipSynThresholds(synThresholdD,leafStrain,tree)
            # since strainPair[0] and strainPair[1] same, no need to reverse
            synThresholdD['minSynThreshold'][strainPair] = minSynThreshold
            synThresholdD['synAdjustThreshold'][strainPair] = synAdjustThreshold
    return synThresholdD
            
def getSynThresholds(binHeightL,indexToBinCenterL,paramD):
    '''Calculate and return minSynThresh and synAdjustThreshold.'''

    binWidth = 1.0/paramD['scoreHistNumBins'] # since scores range from 0-1
    
    # get synteny peak
    peakL = []

    # case 1 (High prominence, narrow peak. Close relatedness
    # in order to get a rightmost peak, if any, we add a dummy bin of
    # height 0 on right
    tempBinHeightL = numpy.append(binHeightL,0)
    tempIndexToBinCenterL = numpy.append(indexToBinCenterL,1)

    L = findPeaksOneCase(tempBinHeightL,tempIndexToBinCenterL,binWidth,paramD['synPeakWidthCase1'],paramD['widthRelHeight'],paramD['synRequiredProminenceCase1'],paramD['synLeftPeakLimit'],paramD['synRightPeakLimit'])
    peakL.extend(L)
    
    # case 2 (wide width, low prominence. Distant relatedness.)
    L = findPeaksOneCase(binHeightL,indexToBinCenterL,binWidth,paramD['synPeakWidthCase2'],paramD['widthRelHeight'],paramD['synRequiredProminenceCase2'],paramD['synLeftPeakLimit'],paramD['synRightPeakLimit'])
    peakL.extend(L)
    
    if peakL == []:
        minSynThreshold = paramD['defaultMinSynThreshold']
        synAdjustThreshold = paramD['defaultSynAdjustThreshold']
    else:
        # sort on position, in case more than one peak (unlikely), take leftmost.
        peakL.sort(key=lambda x: x[1])

        peakHeight,peakPos,leftExtremeOfPeakPos,rightExtremeOfPeakPos = peakL[0]
        minSynThreshold = leftExtremeOfPeakPos - paramD['minSynThresholdIncrement']
        synAdjustThreshold = leftExtremeOfPeakPos - paramD['synAdjustThresholdIncrement']

    return minSynThreshold,synAdjustThreshold

def getTipSynThresholds(synThresholdD,leafStrain,tree):
    '''Get the synteny threshold values at a tip called leafStrain by
averaging the values between that strain and its nearest neighbors.'''
    
    neighbL = trees.getNearestNeighborL(leafStrain,tree)
    avSyn = 0
    avSynAdjust = 0
    for neighb in neighbL:
        avSyn += synThresholdD['minSynThreshold'][neighb,leafStrain]
        avSynAdjust += synThresholdD['synAdjustThreshold'][neighb,leafStrain]
    avSyn /= len(neighbL)
    avSynAdjust /= len(neighbL)
    return avSyn,avSynAdjust


#### Family creation functions

def createAllFamiliesDescendingFromInternalNode(subtreeL,familyMrca,nodeGenesL,scoresO,geneNames,absMinRawThresholdForHomologyD,synThresholdD,paramD,geneUsedL,familiesO,famNumCounter,locusFamNumCounter):
    '''Creates all Families and subsidiary LocusFamilies descending from
the node rooted familyMrca. Basic parts of the Phigs algorithm are
here. Creating the seeds, and using them to get a family. (With very
minor changes in the use of the norm score and synteny
adjustment). But then we divide the family into LocusFamilies, which
is not Phigs.
    '''

    subtree=subtreeL[familyMrca]
    
    # for use in createAllLocusFamiliesDescendingFromInternalNode call below
    familySubtreeNodeOrderL = createNodeProcessOrderList(subtree)

    leftS,rightS = createLRSets(subtreeL,familyMrca,nodeGenesL,None)

    seedL = createSeedL(leftS,rightS,scoresO,geneNames,absMinRawThresholdForHomologyD,paramD)
    for seed in seedL:
        # each seed corresponds to a prospective gene family.
        seedRawSc,seedG1,seedG2 = seed

        if seedRawSc == -float('inf'):
            # we've gotten to the point in the seed list with
            # genes having no match on the other branch
            break
        else:
            # getting initial family, using only raw score and synteny bump
            famS=createFamilyFromSeed(seedG1,seedG2,geneUsedL,scoresO,leftS,rightS,geneNames,seedRawSc,absMinRawThresholdForHomologyD,synThresholdD,paramD)

            if famS == None:
                # one of the genes the the family was already
                # used, so createFamilyFromSeed returned None
                continue
            else:
                # none of the genes in famS used yet
                for gene in famS:
                    geneUsedL[gene] = True

                # now set up familiesO to take this family and
                # determine the corresponding locusFamilies
                familiesO.initializeFamily(famNumCounter,familyMrca,[seedG1,seedG2])
                locusFamNumCounter,familiesO = createAllLocusFamiliesDescendingFromInternalNode(subtreeL,familyMrca,geneNames,famS,[seedG1,seedG2],famNumCounter,locusFamNumCounter,scoresO,paramD,synThresholdD,familiesO,familySubtreeNodeOrderL,nodeGenesL)
                famNumCounter+=1 # important to increment this after call to createAllLocusFamiliesDescendingFromInternalNode
                
    return geneUsedL,locusFamNumCounter,famNumCounter,familiesO

def createAllFamiliesAtTip(nodeGenesL,familyMrca,geneUsedL,tipFamilyRawThresholdD,scoresO,geneNames,absMinRawThresholdForHomologyD,synThresholdD,paramD,familiesO,famNumCounter,locusFamNumCounter):
    '''Creates all Families and subsidiary LocusFamilies at the tip
given by familyMrca. Because we've come through the nodes in
pre-order, we know that all unused genes at this node are in a
families with mrca here. (they can still be multi gene families).
    '''

    unusedGenesAtThisTipS=set()
    for gene in nodeGenesL[familyMrca]: # familyMrca is a tip
        # gene is at this tip
        if not geneUsedL[gene]:
            # not used yet
            unusedGenesAtThisTipS.add(gene)

    # pull out the threshold we'll use given the strain
    tipFamilyRawThreshold = tipFamilyRawThresholdD[familyMrca]

    # Now we pull out families greedily. Simply take a gene,
    # and get all other genes that isSameFamily says are shared
    while len(unusedGenesAtThisTipS)>0:
        seed = unusedGenesAtThisTipS.pop()
        newFamS=set([seed])
        for newGene in unusedGenesAtThisTipS:

            if scoresO.isEdgePresentByEndNodes(seed,newGene):

                addIt = isSameFamily(seed,newGene,scoresO,geneNames,tipFamilyRawThreshold,absMinRawThresholdForHomologyD,synThresholdD,paramD)
                
                if addIt:
                    newFamS.add(newGene)

        # newFamS now contains a set of genes with significant
        # similarity. Remove it from unusedGenesAtThisTipS,
        # and create a new family from it.
        unusedGenesAtThisTipS.difference_update(newFamS)
        for gene in newFamS: # mark these as used
            geneUsedL[gene] = True


        familiesO.initializeFamily(famNumCounter,familyMrca)
        lfOL,locusFamNumCounter = createAllLocusFamiliesAtTip(newFamS,geneNames,familyMrca,scoresO,paramD,synThresholdD,famNumCounter,locusFamNumCounter)
        
        for lfO in lfOL:
            familiesO.addLocusFamily(lfO)

        famNumCounter+=1 # important to increment this after creating LocusFamilies

    return geneUsedL,locusFamNumCounter,famNumCounter,familiesO

def createLRSets(subtreeL,mrca,nodeGenesL,restrictS):
    '''At given mrca, obtain all genes in species in left branch and put
in leftS, and all genes from species in right branch to
rightS. Restrict each of these to be only genes in restrictS. If
restrictS is None, then use all genes.
    '''

    subtree=subtreeL[mrca]
    
    leftS=set()
    for tip in trees.leafList(subtree[1]):
        leftS.update(nodeGenesL[tip])
    rightS=set()
    for tip in trees.leafList(subtree[2]):
        rightS.update(nodeGenesL[tip])

    if restrictS != None:    
        leftS.intersection_update(restrictS)
        rightS.intersection_update(restrictS)
        
    return(leftS,rightS)

def closestMatch(gene,S,scoresO,geneNames,absMinRawThresholdForHomologyD,paramD):
    '''Find the closest match to gene among the genes in the set S in the
graph scoresO. Eliminate any matches that have a raw score below what
is in homologyHomologyRawThresholdD, a coreSynSc below
minCoreSynThresh, or a synteny score below synThresholdD.

    '''
    bestGene=None
    bestEdgeScore = -float('inf')
    for otherGene in scoresO.getConnectionsGene(gene):
        if otherGene in S:
            if isSameFamily(gene,otherGene,scoresO,geneNames,bestEdgeScore,absMinRawThresholdForHomologyD,None,paramD):
                # we don't want to use synThresholdD, hence the Nones
                bestEdgeScore = scoresO.getScoreByEndNodes(gene,otherGene,'rawSc')
                bestGene = otherGene
    return bestEdgeScore, gene, bestGene
    
def createSeedL(leftS,rightS,scoresO,geneNames,absMinRawThresholdForHomologyD,paramD):
    '''Create a list which has the closest match for each gene on the
opposite side of the tree. e.g. if a gene is in tree[1] then we're
looking for the gene in tree[2] with the closest match. We eliminate
any matches that are below threshold for normalized or syntenty
scores. For each gene we get this closest match, put in a list, sort,
and return.
    '''
    seedL=[]
    for gene in leftS:
        seedL.append(closestMatch(gene,rightS,scoresO,geneNames,absMinRawThresholdForHomologyD,paramD))
    for gene in rightS:
        seedL.append(closestMatch(gene,leftS,scoresO,geneNames,absMinRawThresholdForHomologyD,paramD))
    seedL.sort(reverse=True)
    return seedL

def createFamilyFromSeed(g1,g2,geneUsedL,scoresO,leftS,rightS,geneNames,thisFamRawThresh,absMinRawThresholdForHomologyD,synThresholdD,paramD):
    '''Based on a seed (seedScore, g1, g2) search for a family. Using the
PhiGs approach, we collect all genes which are closer to members of
the family than the two seeds are from each other. We have a normScore
threshold below which we will not add genes. We also have a synteny
adjustment of the raw score where we make the raw score between a pair
a bit better if their synteny is above the species pair specific
adjustThresh in synThresholdD. In general, if a gene has syntenic
connections to genes already in the family, this makes us more
confident that this gene belongs in the family. Returns a set
containing genes in the family.

    '''

    if geneUsedL[g1] or geneUsedL[g2]:
        # one of these has been used already, stop now.
        return None

    famS = set()
    genesToSearchForConnectionsS = set([g1,g2])

    while len(genesToSearchForConnectionsS) > 0:
        
        matchesS = getFamilyMatches(genesToSearchForConnectionsS,scoresO,leftS,rightS,famS,geneNames,thisFamRawThresh,absMinRawThresholdForHomologyD,synThresholdD,paramD,geneUsedL)

        if matchesS == None:
            return None
        
        famS.update(genesToSearchForConnectionsS)
        genesToSearchForConnectionsS = matchesS
        
    return famS

def getFamilyMatches(genesToSearchForConnectionsS,scoresO,leftS,rightS,famS,geneNames,thisFamRawThresh,absMinRawThresholdForHomologyD,synThresholdD,paramD,geneUsedL):
    ''''''
    matchesS=set()
    for famGene in genesToSearchForConnectionsS:
        for newGene in scoresO.getConnectionsGene(famGene):
            if newGene in leftS or newGene in rightS:
                # it is from a species descended from the node
                # we're working on
                if newGene not in genesToSearchForConnectionsS and newGene not in famS:
                    # it shouldn't been in our current list to search,
                    # or be one we've already put in the family (or
                    # we'll waste effort)
                    addIt = isSameFamily(famGene,newGene,scoresO,geneNames,thisFamRawThresh,absMinRawThresholdForHomologyD,synThresholdD,paramD)
                    
                    if addIt:
                        if geneUsedL[newGene]:
                            # this one's been used already. That
                            # means the whole family should be
                            # thrown out. Just stop now.
                            return None
                        else:
                            matchesS.add(newGene)
    return matchesS

def isSameFamily(famGene,newGene,scoresO,geneNames,thisFamRawThresh,absMinRawThresholdForHomologyD,synThresholdD,paramD):
    '''Given famGene that is inside a family, and newGene we are
considering adding, check the various scores to determine if we should
add it. Return boolean.
    '''

    rawSc = scoresO.getScoreByEndNodes(famGene,newGene,'rawSc')
    synSc = scoresO.getScoreByEndNodes(famGene,newGene,'synSc')

    # get minThresh from absMinRawThresholdForHomologyD
    strain1 = geneNames.numToStrainNum(famGene)
    strain2 = geneNames.numToStrainNum(newGene)
    absoluteMinRawThresh = absMinRawThresholdForHomologyD[(strain1,strain2)]

    if synThresholdD == None:
        synAdjustThresh = float('inf')
    else:
        # change later
        synAdjustThresh = synThresholdD['synAdjustThreshold'][(strain1,strain2)]

        
    addIt = False
    if rawSc < absoluteMinRawThresh:
        # Raw score is simply too low, don't add newGene. Raw score is
        # below the minimum value we've calculated for this species
        # pair.  (Modification of PhiGs)
        pass
    elif rawSc >= thisFamRawThresh:
        # If its within the seed distance, add it
        # (basic PhiGs approach). we have the =
        # there in case thisFamRawThresh is 1.
        addIt = True
    elif synSc > synAdjustThresh:
        # its above the syn score adjustment
        # threshold, so increase rawSc a bit. This
        # addresses a problem with closely related
        # families where the seed score is very
        # similar. Sometimes by chance things
        # that should have been added weren't
        # because they weren't more similar than
        # an already very similar seed.
        # (Modification of PhiGs)
        adjSc = rawSc * paramD['synAdjustExtent']
        if adjSc > 1: adjSc = 1 # truncate back to 1
        if adjSc >= thisFamRawThresh:
            addIt = True
                
    return addIt

def createAllLocusFamiliesDescendingFromInternalNode(subtreeL,familyMrca,geneNames,famGenesToSearchS,seedPairL,famNumCounter,locusFamNumCounter,scoresO,paramD,synThresholdD,familiesO,familySubtreeNodeOrderL,nodeGenesL):
    '''Given a family in famGenesToSearchS, break it up into subsidiary locus families
based on synteny. We iterate through the subtree rooted at familyMrca
in pre-order (ancestors first). Using seeds, we try to find groups
among famS that share high synteny.'''

    # split out LocusFamilies at non-syntenic locations
    for lfMrca,lchild,rchild in familySubtreeNodeOrderL:

        if lchild != None:
            # not a tip
            lfOL,locusFamNumCounter,famGenesToSearchS = createAllLocusFamiliesAtOneInternalNode(subtreeL,lfMrca,nodeGenesL,geneNames,famGenesToSearchS,scoresO,paramD,synThresholdD,famNumCounter,locusFamNumCounter)
            
            for lfO in lfOL:
                familiesO.addLocusFamily(lfO)

        else:
            # we're at a tip.
            
            # get only the genes at this tip
            genesAtThisTipS = famGenesToSearchS.intersection(nodeGenesL[lfMrca])

            # remove them from famGenesToSearchS
            famGenesToSearchS.difference_update(genesAtThisTipS)

            # Get lf objects for all these genes
            lfOL,locusFamNumCounter = createAllLocusFamiliesAtTip(genesAtThisTipS,geneNames,lfMrca,scoresO,paramD,synThresholdD,famNumCounter,locusFamNumCounter)            

            # add to our families object
            for lfO in lfOL:
                familiesO.addLocusFamily(lfO)

    return locusFamNumCounter,familiesO

def createAllLocusFamiliesAtOneInternalNode(subtreeL,lfMrca,nodeGenesL,geneNames,famGenesToSearchS,scoresO,paramD,synThresholdD,famNumCounter,locusFamNumCounter):
    '''Obtains all locus families at the internal node defined by lfMrca.'''

    lfOL = []
    while True:

        lfSeedPairL = createLFSeed(subtreeL,lfMrca,nodeGenesL,geneNames,famGenesToSearchS,scoresO,paramD,synThresholdD)
        
        if lfSeedPairL == []:
            # there are no (more) seeds stradling this internal node,
            # break out
            break

        lfO,locusFamNumCounter,famGenesToSearchS = createLocusFamilyFromSeed(famNumCounter,locusFamNumCounter,lfMrca,lfSeedPairL,famGenesToSearchS,subtreeL,geneNames,scoresO,paramD,synThresholdD)
        lfOL.append(lfO)
        
    return lfOL,locusFamNumCounter,famGenesToSearchS
            
def createLFSeed(subtreeL,lfMrca,nodeGenesL,geneNames,famGenesToSearchS,scoresO,paramD,synThresholdD):
    '''Given a set of genes famGenesToSearchS from a family, try to find a
seed based at lfMrca. A seed consists of two genes, one in the left
subtree and one in the right, which are syntenically consistent.
    '''
    
    leftS,rightS = createLRSets(subtreeL,lfMrca,nodeGenesL,famGenesToSearchS)
    
    for lGene in leftS:
        for rGene in rightS:
            if isSameLocusFamily(lGene,rGene,scoresO,geneNames,paramD,synThresholdD):
                return [lGene,rGene]
    return []
        
def createLocusFamilyFromSeed(famNumCounter,locusFamNumCounter,lfMrca,seedPairL,famGenesToSearchS,subtreeL,geneNames,scoresO,paramD,synThresholdD):
    '''Returns a LocusFamily object, containing genes associated with
those in seedPairL, in the subtree definied at lfMrca. Does single
linkage clustering, adding in anything in famGenesToSearchS with above
threshold synteny. Note that these seeds are not the seed from family formation (which might not be syntenic) but rather an independently generated pair which we know belong in the same LocusFamily
    '''

    lfO = LocusFamily(famNumCounter,locusFamNumCounter,lfMrca)
    locusFamNumCounter+=1
    
    famGenesToSearchS.difference_update(seedPairL)
    lfO.addGenes(seedPairL,geneNames)

    subtree=subtreeL[lfMrca]
    strainL = trees.leafList(subtree)
        
    while True:
        genesToAddS = getLocusFamilyMatches(lfO,famGenesToSearchS,geneNames,strainL,scoresO,paramD,synThresholdD)
        if len(genesToAddS) == 0:
            break

        famGenesToSearchS.difference_update(genesToAddS)
        lfO.addGenes(genesToAddS,geneNames)
        
    return lfO,locusFamNumCounter,famGenesToSearchS

def getLocusFamilyMatches(lfO,famGenesToSearchS,geneNames,strainL,scoresO,paramD,synThresholdD):
    '''Given a LocusFamily object lfO and some remaining genes, search
through the remaining genes to find those that match syntenically and
are in a child species of lfMrca. Return a list of genes to add.
    '''
    genesToAddS=set()
    for searchGene in famGenesToSearchS:
        # test if searchGene is in a child species of lfMrca
        if geneNames.numToStrainNum(searchGene) in strainL:
            # this searchGene is in a strain that is a child of the lfMrca we're working on
            for lfGene in lfO.iterGenes():
                # we don't use absMinRawThreshold, thisFamRawThresh, or
                # synAdjustThresh. If the pair have values above
                # minCoreSynThresh and minSynThres, then addIt will be
                # True.
                addIt = isSameLocusFamily(searchGene,lfGene,scoresO,geneNames,paramD,synThresholdD)
                if addIt:
                    genesToAddS.add(searchGene)
                    break
                
    return genesToAddS

def isSameLocusFamily(gene1,gene2,scoresO,geneNames,paramD,synThresholdD):
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

    strain1 = geneNames.numToStrainNum(gene1)
    strain2 = geneNames.numToStrainNum(gene2)
    minSynThresh = synThresholdD['minSynThreshold'][(strain1,strain2)]

    if coreSynSc < paramD['minCoreSynThresh'] or synSc < minSynThresh:
        # one of the two types of synteny below threshold, so this
        # pair doesn't meet the requirements for being in the same
        # LocusFamily
        addIt = False
    else:
        addIt = True
            
    return addIt

def createAllLocusFamiliesAtTip(genesAtThisTipS,geneNames,lfMrca,scoresO,paramD,synThresholdD,famNumCounter,locusFamNumCounter):
    '''Given a set of genes famGenesToSearchS, search for all those found
on the tip given by lfMrca. Break these into LocusFamilies. Many will
be single gene LocusFamilies, but some may be multi-gene
    '''

    # Break these up into LocusFamilies. Many will be
    # single gene LocusFamilies, but some may be multi-gene
    lfGroupsL=[]
    while len(genesAtThisTipS) > 0:

        seed = genesAtThisTipS.pop()

        currentGroupS=set([seed])
        for gene in genesAtThisTipS:
            addIt = isSameLocusFamily(seed,gene,scoresO,geneNames,paramD,synThresholdD)
            if addIt:
                currentGroupS.add(gene)

        genesAtThisTipS.difference_update(currentGroupS)

        lfGroupsL.append(currentGroupS)

    lfOL=[]
    for lfGroupS in lfGroupsL:

        lfO = LocusFamily(famNumCounter,locusFamNumCounter,lfMrca)
        locusFamNumCounter+=1
        lfO.addGenes(lfGroupS,geneNames)
        lfOL.append(lfO)

    return lfOL,locusFamNumCounter

## Input/output

def writeFamilies(familiesO,geneNames,strainNum2StrD,fileName):
    '''Write all gene families to fileName, one family per line.'''
    f=open(fileName,'w')
    for fam in familiesO.iterFamilies():
        f.write(fam.fileStr(strainNum2StrD,geneNames)+'\n')
    f.close()


def readFamilies(familyFN,tree,geneNames,strainStr2NumD):
    '''Read the family file named familyFN, creating a Families object.
    '''

    familiesO = Families(tree)
    f=open(familyFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.split('\t')
        famNum=int(L[0])
        mrca = strainStr2NumD[L[1]]
        if L[2] == "-":
            seedPairL = None
        else:
            seedPairL = [L[2],L[3]]

        lfL = L[4:]

        familiesO.initializeFamily(famNum,mrca,seedPairL)

        for lfStr in lfL:
            lfSplitL = lfStr.rstrip().split(',')
            locusFamNum=int(lfSplitL[0])
            lfMrca = strainStr2NumD[lfSplitL[1]]
            geneL=[]
            for geneName in lfSplitL[2:]:
                geneL.append(geneNames.nameToNum(geneName))
            lfO = LocusFamily(famNum,locusFamNum,lfMrca)
            lfO.addGenes(geneL,geneNames)
            familiesO.addLocusFamily(lfO)

    f.close()
    return familiesO
