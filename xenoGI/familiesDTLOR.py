import sys,numpy,os,random,glob,time
from scipy.signal import find_peaks
from Bio import Phylo
from multiprocessing import Pool
from . import trees,scores,DTLOR_DP,treeParser # consolidate treeParser with trees later
from .Family import *
from .analysis import printTable
#import matplotlib.pyplot as plt
from collections import deque

#### Main function

def createFamiliesO(tree,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,outputSummaryF):
    '''Main function to create families. First creates an initial families
object, then reconciles those against the species tree to make an
originFamilies object.
    '''

    # threshold for maximum size of gene tree, if above, split
    maxFamilySize=60

    # define some variables
    initFamilyFN = paramD['initFamilyFN']
    originFamilyFN =  paramD['originFamilyFN']
    geneInfoFN = paramD['geneInfoFN']

    # checks
    homologyCheck(genesO,aabrhHardCoreL,scoresO,outputSummaryF,paramD)
    
    # create initial families
    initialFamiliesO, locusMapD=createInitialFamiliesO(maxFamilySize,tree,scoresO,genesO,aabrhHardCoreL,paramD)
    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT)) # needed for writeFamilies
    writeFamilies(initialFamiliesO,initFamilyFN,genesO,strainNamesT,paramD)
    print("Initial families:",file=outputSummaryF)
    writeFamilyFormationSummary(initialFamiliesO,outputSummaryF)
    
    print("Initial families written out to %s"%initFamilyFN,file=sys.stderr)

    # make gene family trees
    trees.makeGeneFamilyTrees(paramD,genesO,initialFamiliesO) #create gene tree for each initial family 
    print("Finished making gene trees")

    #initialFamiliesO = readFamilies(paramD['initFamilyFN'],tree,genesO) # TEMP

    # create origin families
    originFamiliesO = reconcile(tree,initialFamiliesO,locusMapD,genesO,paramD,strainNamesT)
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

def createInitialFamiliesO(maxFamilySize,tree,scoresO,genesO,aabrhHardCoreL,paramD):
    '''Given a scoresO object, create initial genes families based on
    blast connectivity and return as a familiesO object.
    '''
    
    # get synteny thesholds for locus family formation
    synThresholdD = getSynThresholdD(paramD,scoresO,genesO,aabrhHardCoreL,tree)

    # get initial families as list of sets
    initFamilySetL, degreeD = createInitialFamilySetL(scoresO, genesO)
    print("Number of initial families before size control:",len(initFamilySetL),file=sys.stderr)
    
    initFamilySetL.sort(reverse=True,key=len)  #sort by number of genes in descending order
    oversizedFamiliesL = [familyS for familyS in initFamilySetL if len(familyS)>maxFamilySize]
    print("Number of initial families with more than %d genes: %d"%(maxFamilySize,len(oversizedFamiliesL)))
    # for index, family in enumerate(oversized): 
    #     visualize_connectivity(family,"oversized_family_%d"%index,scoresO,degreeD)

    # divide these families into locus families, and split oversized
    initFamilyDivIntoLocFamsL=[]
    for i,familyS in enumerate(oversizedFamiliesL):
        # visualize_connectivity(familyS,"oversizedFam_%d"%i,scoresO, degreeD)
        locusFamLL=divideInitialFamilyIntoLocusFamilies(familyS, genesO,scoresO,paramD,synThresholdD)  #list of lists
        splitFamsL=splitOversizedFamByLocus(familyS, scoresO,maxFamilySize, locusFamLL)
        for j, (familyS, locusFamLL) in enumerate(splitFamsL):
            # visualize_connectivity(familyS,"%i_postSplit_%d"%(i,j),scoresO)
            initFamilyDivIntoLocFamsL.append((len(familyS),familyS,locusFamLL))
            
    for familyS in initFamilySetL[len(oversizedFamiliesL):]:
        locusFamLL=divideInitialFamilyIntoLocusFamilies(familyS, genesO,scoresO,paramD,synThresholdD)
        initFamilyDivIntoLocFamsL.append((len(familyS),familyS,locusFamLL))
      
    # sort the list, descending by size (so treeFNs in dir are thus ordered)
    initFamilyDivIntoLocFamsL.sort(reverse=True,key=lambda x: x[0])

    # create output objects
    locusMapD={}
    # construct Family object for initial families
    initFamiliesO = Families(tree)

    # store in familiesO object
    # these locusFam Ids need to be unique across all, start counting from 1. EB: WHY not 0?
    famNumCounter=1
    locusFamNumCounter=1
    totalAddedTolocusFamilies=0
    for size, familyS, locusFamLL in initFamilyDivIntoLocFamsL:
        locusFamNumCounter, famNumCounter,totalAddedTolocusFamilies=addFamilyToInitialFamiliesO(familyS,initFamiliesO,famNumCounter,locusFamLL,locusMapD,locusFamNumCounter,genesO,tree,totalAddedTolocusFamilies)
    return initFamiliesO, locusMapD

def getSynThresholdD(paramD,scoresO,genesO,aabrhHardCoreL,tree):
    '''Creates a dictionary to store synteny thresholds. This dictionary
itself contains two dictionaries one each for minCoreSyntThresh
(minimum core synteny score allowed for family formation), and
minSynThresh (minimum synteny score allowed for family
formation). These dictionaries in turn are keyed by strain pair.
    '''
    quantileForObtainingSynThresholds = paramD['quantileForObtainingSynThresholds']
    multiplierForObtainingSynThresholds = paramD['multiplierForObtainingSynThresholds']
    
    synThresholdD = {}
    synThresholdD['minSynThreshold'] = {}
    synThresholdD['minCoreSynThreshold'] = {}

    # coreSynSc
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
            synThresholdD['minCoreSynThreshold'][strainPair] = getTipThreshold(tree,leafStrain,synThresholdD,'minCoreSynThreshold')
            synThresholdD['minSynThreshold'][strainPair] = getTipThreshold(tree,leafStrain,synThresholdD,'minSynThreshold')

    return synThresholdD

def getTipThreshold(tree,leafStrain,synThresholdD,thresholdType):
    '''Get the synteny threshold values at a tip called leafStrain by
averaging the values between that strain and its nearest neighbors.'''
    
    neighbL = trees.getNearestNeighborL(leafStrain,tree)
    avThresh = 0
    for neighb in neighbL:
        strainPair = tuple(sorted((neighb,leafStrain)))
        avThresh += synThresholdD[thresholdType][strainPair]
    avThresh /= len(neighbL)
    return avThresh

def createInitialFamilySetL(scoresO, genesO):
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
    degreeD={}
   
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
    allGenesL=list(genesO.iterGenes())
    scoresO.createNodeConnectD() 
    connectedGenes=list(scoresO.nodeConnectD.keys())
    initFamilySetL=[]
    visited=set()
    for gene in allGenesL: 
        if gene in scoresO.nodeConnectD:
            if gene not in visited: 
                temp =set()
                newFam=stronglyConnected(temp, gene, visited)
                
                initFamilySetL.append(newFam) 
        else:
            fam=set()
            fam.add(gene)
            initFamilySetL.append(fam)
  
    return initFamilySetL, degreeD

def addFamilyToInitialFamiliesO(familyS,initFamiliesO,famNumCounter,locusFamLL,locusMapD,locusFamNumCounter,genesO,tree,totalAddedTolocusFamilies):
    '''Given a families object, a family set and a list of locus families
in that family, add the family to the object.
    '''
    speciesL=list(set([genesO.numToStrainName(gene) for gene in familyS]))
    mrca=findMrca(speciesL, tree)
    #add each initial family as a Family object (still empty)
    initFamiliesO.initializeFamily(famNumCounter,mrca) 
    for locusFamilyL in locusFamLL: #locusFamilyL contains all the genes in that lf
        speciesL=[]
        for gene in locusFamilyL:
            locusMapD[str(gene)]=locusFamNumCounter
            speciesL.append(genesO.numToStrainName(gene))
        lfMrca=findMrca(speciesL, tree)
        lf=LocusFamily(famNumCounter,locusFamNumCounter,lfMrca)
        lf.addGenes(locusFamilyL, genesO)
        totalAddedTolocusFamilies+=len(locusFamilyL)
        initFamiliesO.addLocusFamily(lf)
        locusFamNumCounter+=1
    famNumCounter+=1
    return locusFamNumCounter, famNumCounter,totalAddedTolocusFamilies

def splitOversizedFamByLocus(family, scoresO,maxFamilySize, locusFamLL):
    """Split families that are too large without breaking up locus
    families.  Returns a list of new smaller families with their locus
    families in (family, [locus families]) pair
    """
    origin_size=len(family)
    locusFamLL.sort(reverse=True,key=len)
    #calculate pairwise RawSc between locus families
    similarity={}
    for i in range(len(locusFamLL)):
        loc1=locusFamLL[i]
        if loc1[0] not in similarity:
            similarity[loc1[0]]={}
        for j in range(i+1,len(locusFamLL)):
            loc2=locusFamLL[j]
            if loc2[0] not in similarity:
                similarity[loc2[0]]={}
            maxSc=0
            for gene1 in loc1:
                for gene2 in loc2:
                    if scoresO.isEdgePresentByEndNodes(gene1, gene2):
                        rawSc=scoresO.getScoreByEndNodes(gene1, gene2, "rawSc")
                    else: rawSc=0
                    if rawSc>maxSc: maxSc=rawSc
            similarity[loc1[0]][loc2[0]]=maxSc  #use the first gene in the loc as identifier
            similarity[loc2[0]][loc1[0]]=maxSc
    #fill the buckets greedily
    allLfBucketsL=[]
    current_bucket=[]
    currentSize=0
    current_loc=None
    while locusFamLL!=[]: #while there are still more locus to add
        if current_loc==None:
            current_loc=locusFamLL[0]  #start with the biggest loc fam in queue
        
        else: #find the loc with the highest similarity score and not in any bucket yet
            #current_loc should already be deleted form the list
            max_ind=numpy.argmax([similarity[current_loc[0]][other_loc[0]] for other_loc in locusFamLL])
            current_loc=locusFamLL[max_ind]
        if len(current_loc)>maxFamilySize:
            print("There is a locus family with size over threshold, putting that as one initial family")
        #try putting into the current bucket
        if (currentSize==0) or (currentSize+len(current_loc)<=maxFamilySize):
            #always allow putting into empty bucket, or if adding to bucket doesn't overflow
            current_bucket.append(current_loc)
            currentSize+=len(current_loc)
            locusFamLL.remove(current_loc)
        else:#start a new bucket
            allLfBucketsL.append(current_bucket) #archive the current one
            current_bucket=[]
            currentSize=0
            current_loc=None
    allLfBucketsL.append(current_bucket) #archive the current one

    splitFamsL=[]
    newSize=0
    for lfBucketL in allLfBucketsL: # lfBucketL is a list of locus families
        newInitialFamilyS=set((gene for locus_family in lfBucketL for gene in locus_family))
        newSize+=len(newInitialFamilyS)
        splitFamsL.append((newInitialFamilyS,lfBucketL))

    return splitFamsL

def visualize_connectivity(initial_family,fn,scoresO,degreeD=None):
    if degreeD==None: #not the same as the whole gene set 
        degreeD={}
        genes=list(initial_family)
        for gene in genes:
            all_neighbors=set(scoresO.getConnectionsGene(gene))
            neighbors_in_family=[neighbor for neighbor in genes if neighbor in all_neighbors]
            degreeD[gene]=len(neighbors_in_family)
    degrees=[degreeD[gene] for gene in initial_family]
    plt.hist(degrees,bins=20, color='b', edgecolor='k', alpha=0.65)
    plt.axvline(sum(degrees)/float(len(degrees)), color='k', linestyle='dashed', linewidth=1, label='Mean: {:.2f}'.format(sum(degrees)/float(len(degrees))))
    plt.axvline(len(degrees)/float(2), color='c', linestyle='dotted', linewidth=1, label='N/2: {:.2f}'.format(len(degrees)/float(2)))
    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    # min_ylim, max_ylim = plt.ylim()
    # plt.text(sum(degrees)/float(len(degrees))*1.01, max_ylim*0.7, 'Mean: {:.2f}'.format(sum(degrees)/float(len(degrees))))   
    # plt.text(len(degrees)/float(2)*1.01, max_ylim*0.9, 'Half of nodes: {:.2f}'.format(len(degrees)/float(2))) 
    plt.gca().set(title='%s Histogram'%fn, ylabel='frequency', xlabel="Number of neighbors by BLAST")
    if not os.path.exists("initialFamily_histogram"):
        os.mkdir("initialFamily_histogram")

    plt.savefig("initialFamily_histogram/%s"%fn)
    plt.close()


#THIS CURRENTLY RUNS VERY SLOW DUE TO THE STRUCTURE OF THE TREE
def findMrca(speciesL, tree):
    '''Given a list of species and a tree, return the node of the most recent common ancestor.'''
    
    #find the LCA recursively 
    #base cases

    if len(speciesL)==1:
        return speciesL[0]

    species1=speciesL[0]
    species2=speciesL[1]

    if len(speciesL)==2:
        return findMrcaPair(species1,species2,tree)
    
    else:
        LCA=findMrcaPair(species1, species2, tree)

        newSpecies=[LCA]
        newSpecies.extend(speciesL[2:])

        return findMrca(newSpecies, tree)

def findMrcaPair(species1, species2, tree):
    '''Return node of most recent common ancestor of two species.'''
    if len(tree)==0:
        return None
    elif tree[0]==species1 or tree[0]==species2:
        #if root is equal to either, return root as the lCA
        return tree[0]
    else:
        left=findMrcaPair(species1,species2, tree[1])
        right=findMrcaPair(species1, species2, tree[2])
        if (left is not None) and (right is not None):
            #they are in different subtrees, return the root
            return tree[0]
        if left is None:
            return right
        else: return left

def divideInitialFamilyIntoLocusFamilies(familyS, genesO, scoresO, paramD,synThresholdD):
    """
    Input 
    ----------------------------------------
    familyS:            An initial family of genes generated by 'createInitialFamilySetL'. This is represented 
                        by a set that contains genes that are connected by BLAST results
    scoresO:            The Score object that contains the information for all pairwise synteny scores

    Return 
    ----------------------------------------
    locusFamilies:      A list of lists, sub-lists contain locus families produced by clustering

    """

    minGenes = 2 # don't divide if num genes less than this.

    
    #if there are too few genes, we should probably not cluster further
    if len(familyS)<=minGenes:
        locusFamilyL=list(familyS)
        return [locusFamilyL]
    #construct the fully connected graph with weighted edges 
    num_gene=len(familyS)
    genesAr=numpy.array(list(familyS))
    geneNeighborsD={}
    #note there is no edge going from one to self.
    
    #loop over all pairs of genes
    for i in range(len(genesAr)):
        for j in range(i+1, len(genesAr)):
            gene1=genesAr[i]
            gene2=genesAr[j]
            if scoresO.isEdgePresentByEndNodes(gene1,gene2):
                sameLocusFam = isSameLocusFamily(gene1,gene2,scoresO,genesO,paramD,synThresholdD)
                #the graph is symmetrical

                # initialize the adjacency dict if key doesn't already exist
                if gene1 not in geneNeighborsD:
                    geneNeighborsD[gene1]=[]
                if gene2 not in geneNeighborsD:
                    geneNeighborsD[gene2]=[]
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
        if neighborsL:
            for neighbGene in neighborsL: 
                if neighbGene not in visitedS:     
                    # Update the list 
                    tempFamL = stronglyConnected(tempFamL, neighbGene, visitedS, geneNeighborsD) 
        return tempFamL 

    locusFamLL=[]
    visitedS=set()
    for gene in genesAr: 
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

def reconcile(speciesTree,initialFamiliesO,locusMapD,genesO,paramD,strainNamesT):
    '''Reconcile gene family trees to the species tree using the DTLOR algorithm.'''

    originFamilyFN=paramD['originFamilyFN']
    geneFamilyTreesDir = paramD['geneFamilyTreesDir']
    numProcesses = paramD['numProcesses']
    D=float(paramD["duplicationCost"])
    T=float(paramD["transferCost"])
    L=float(paramD["lossCost"])
    O=float(paramD["originCost"])
    R=float(paramD["rearrangeCost"])
    print("The parameters used for reconciliation are: D= %.2f, T=%.2f, L=%.2f, O=%.2f, R=%.2f"%(D,T,L,O,R))
    
    gtFileStem = 'fam'

    allGtFilePath = os.path.join(geneFamilyTreesDir,gtFileStem+'*.tre')

    # new families object to record all the origin families
    originFamiliesO = Families(speciesTree)
    # add single-gene initial family as their own origin family
    for init_family in initialFamiliesO.iterFamilies():
        if len(tuple(init_family.iterGenes())) == 1:
            newFamNum=originFamiliesO.getNumFamilies()
            originFamiliesO.initializeFamily(newFamNum,init_family.mrca)
            for locusFamily in init_family.getLocusFamilies():
                locusFamily.famNum=newFamNum
                locusFamily.locusFamNum=originFamiliesO.getNumLocusFamilies()
                originFamiliesO.addLocusFamily(locusFamily)
         
    allTreeFN_L=list(sorted(glob.glob(allGtFilePath)))
    binaryTreeFN_L=[]
    no=0
    yes=0
    for treeFN in allTreeFN_L:
        bpTree = Phylo.read(treeFN, 'newick', rooted=False)
        try:
            bpTree.root_at_midpoint()  #root arbitrarily for further rerooting 
        except:
            bpTree.root_with_outgroup(bpTree.get_terminals()[0])
        treeParser.tabulate_names(bpTree)   #name internal nodes
        if not treeParser.is_binary(bpTree): 
            #get fam num
            prefix=geneFamilyTreesDir+"/"+gtFileStem
            suffix=".tre"
            famNum=treeFN.replace(prefix,'')
            famNum=famNum.replace(suffix,'')
            famNum=int(famNum.lstrip('0'))
            init_fam=initialFamiliesO.getFamily(famNum)
            newFamNum=len(originFamiliesO.familiesD)
            tuple_geneTree=trees.bioPhyloToTupleTree(bpTree)
            originFamiliesO.initializeFamily(newFamNum,init_fam.mrca,tuple_geneTree)
            for locusFamily in init_fam.getLocusFamilies():
                locusFamily.famNum=newFamNum
                locusFamily.locusFamNum=originFamiliesO.getNumLocusFamilies()
                originFamiliesO.addLocusFamily(locusFamily)
            no+=1
        else:
            yes+=1
            binaryTreeFN_L.append(treeFN)
    print("Number of binary trees")
    print(len(binaryTreeFN_L))
    print("Done adding origin families from nonbinary trees")

    # make list of sets of arguments to be passed to p.map. There
    # should be numProcesses sets.
    argumentL = [([],speciesTree,locusMapD,genesO, D, T, L, O, R) for i in range(numProcesses)]
    #distribute all gene tree files into numProcesses separate processes
    for i,treeFN in enumerate(binaryTreeFN_L):
        argumentL[i%numProcesses][0].append(treeFN)
    
    reconLL=[]
    with Pool(processes=numProcesses) as p:
        # store the results to list as they come in
        for reconL in p.imap_unordered(runDTLORGroup, argumentL):  #each process generates a output list for a group of trees
            reconLL.append(reconL)

    print("finished multiprocessing")
    reconLL=reconLL[::-1]  #reversing the order changes which trees are not adding, which means multiprocess is not where things are wrong
    for reconL in reconLL:
        for optGeneRooting, optMPR, treeFN in reconL:
            #update originFamiliesO object
            addOriginFamily(optMPR, optGeneRooting, originFamiliesO, genesO, treeFN) #we know this is called

    print("Number of genes in origin familiesO: %d"%len(originFamiliesO.getAllGenes()))
    writeFamilies(originFamiliesO,originFamilyFN,genesO,strainNamesT,paramD)
    print("Origin families written out to %s"%(originFamilyFN))

    return originFamiliesO

def getTipMapping(geneTree, genesO):
    """
    Fill out the tip mapping (from gene to species) using the binary search function from genomes
    """
    leaves=trees.leafList(geneTree)
    phi={}
    for leaf in leaves:
        geneNum=leaf
        if isinstance(leaf, str):
            geneNum=int(leaf)
        phi[leaf]=genesO.numToStrainName(geneNum)
    return phi


def parseReconForLocus(reconciliation):
    """
    Input
    ---------
    reconciliation:  DTLOR reconciliation dictionary for the whole gene tree
    Return
    ---------
    gl_map:  a dictionary mapping each gene node to its (top, bottom) loci
    given by the DTLOR recon

    mappingNodesWithGene:  a list of keys in the reconciliation that contains the gene mapping
    """
    gl_map=dict()
    mappingNodesWithGene=dict()
    for key in reconciliation.keys():
        gene=key[0]
        locus_t=key[2]
        locus_b=key[3]
        if gene not in gl_map:
            gl_map[gene]=(locus_t,locus_b)
        if gene not in mappingNodesWithGene:
            mappingNodesWithGene[gene]=[key]
        else:
            mappingNodesWithGene[gene].append(key)
    return gl_map, mappingNodesWithGene

def getMRCAforOR(startNode, tree_dict, geneToLocus,MRCA):
    """
    populates the MRCA map from gene leaf to the gene node of the last O or R event
    from it to the root
    """
    locus_t,locus_b=geneToLocus[startNode]
    if locus_b!=locus_t:
        leaves= treeParser.getLeavesInSubtree(startNode, tree_dict)
        for leaf in leaves:
            MRCA[leaf]=startNode
    if startNode in tree_dict:
        left, right=tree_dict[startNode]
        getMRCAforOR(left, tree_dict, geneToLocus, MRCA)
        getMRCAforOR(right, tree_dict, geneToLocus, MRCA)

def getLocusFamiliesInOrigin(startNode, tree_dict, geneToLocus, origin_num, familiesO, genesO, recon):
    """
    Helper function that adds the locus families to originFamiliesO (familiesO). These locus families
    belong in the origin family object that starts with startNode
    """
    LeaftoMRCA={}
    getMRCAforOR(startNode,tree_dict,geneToLocus,LeaftoMRCA)
    MRCAtoLeaves={}
    for gene, MRCA in LeaftoMRCA.items():
        if MRCA in MRCAtoLeaves:
            MRCAtoLeaves[MRCA].append(gene)
        else:
            MRCAtoLeaves[MRCA]=[gene]
    for startNode in MRCAtoLeaves:
        _,locus_b=geneToLocus[startNode]
        _,speciesBr,_=recon[startNode]
        #new locus family with the current gene node's recon mapping species node/branch as the MRCA
        newLocusNum=familiesO.getNumLocusFamilies()
        # lf=LocusFamily(origin_num,int(locus_b),speciesBr)
        lf=LocusFamily(origin_num,newLocusNum,speciesBr, int(locus_b))
        leaves= MRCAtoLeaves[startNode]
        lf.addGenes(leaves, genesO)
        familiesO.addLocusFamily(lf)

def getPartialRecon(reconciliation,mappingNodes):
    """
    Returns a partial dictionary with only key,values pairs
    corresponding to the sub(gene)tree rooted at startNode. 
    Key is a gene node, value is (species branch, locus)
    """
  
    nodeMap={}
    queue=deque(mappingNodes)
    while len(queue)>0:
        currentMapping=queue.popleft()
        value=reconciliation[currentMapping]
        eventType, child1Mapping, child2Mapping=value
        geneBr, speciesBr, tl, bl=currentMapping
        if eventType!='L':  #if it is a loss event the gene node has not been mapped to the bottom 
            nodeMap[geneBr]=(eventType,speciesBr,bl)
        for child in [child1Mapping,child2Mapping]:
            if child!=(None, None, None, None):
                queue.append(child)
    return nodeMap

def addOriginFamily(reconciliation, geneTree,originFamiliesO,genesO, treeFN):
    
    geneToLocus, geneToMappingNodes=parseReconForLocus(reconciliation)
    tree_dict=treeParser.getTreeDictionary(geneTree,{})
   
    #construct an origin family for each Origin event, a locus family for each R event
    for gene in geneToLocus.keys():
        locus_t,locus_b=geneToLocus[gene]
        if locus_t=='*':
            if  locus_b!='*':  #if you use integer conversion, '*' need to be changed
                startNode=gene  
                #parse the reconciliation to only include the part responsible for origin family
                recon=getPartialRecon(reconciliation, geneToMappingNodes[startNode])
                origin_num=originFamiliesO.getNumFamilies()
                _,speciesBr,_=recon[startNode]
                originFamiliesO.initializeFamily(origin_num,speciesBr,geneTree,recon) #use the species for MRCA
                getLocusFamiliesInOrigin(startNode, tree_dict, geneToLocus, origin_num, originFamiliesO, genesO, recon)
   
 
def runDTLORGroup(argT):
    """
    Run reconciliation a group of tree files
    """
    treeFN_list,speciesTree,locusMapD,genesO, D, T, L, O, R= argT
    outputL = []
    for treeFN in treeFN_list:
        optGeneRooting,optMPR = runDTLOR(treeFN,speciesTree,locusMapD,genesO, D, T, L, O, R )
        outputL.append((optGeneRooting,optMPR, treeFN))
    return outputL

def runDTLOR(treeFN,speciesTree, locusMapD,genesO, D, T, L, O, R ):
    bpTree = Phylo.read(treeFN, 'newick', rooted=False)
    try:
        bpTree.root_at_midpoint()  #root arbitrarily for further rerooting 
    except:
        bpTree.root_with_outgroup(bpTree.get_terminals()[0])
    treeParser.tabulate_names(bpTree)   #name internal nodes
    if treeParser.is_binary(bpTree):
        pass
    else: print(treeFN)
    assert treeParser.is_binary(bpTree)
    locus_map = treeParser.rerootingPruning(bpTree, locusMapD)
    # print("gene to locus map")
    # print([(term.name, locus_map[term.name]) for term in bpTree.get_terminals()])

    # for term in bpTree.get_terminals(): term.name=term.name+"_"+str(locus_map[term.name])  #rename the tips adding the syntenic loc
    # Phylo.draw_ascii(bpTree)
    tuple_geneTree=trees.bioPhyloToTupleTree(bpTree)
    phi=getTipMapping(tuple_geneTree,genesO)  #gene to species mapping for these specific gene tree
    
    all_rootings=treeParser.get_all_rerootings(tuple_geneTree, locus_map)
    if all_rootings==[]:  #all rerooting not valid (all nodes have the same loc)
        all_rootings=[tuple_geneTree]
    best_score=float('inf')
    bestMPRs=[]
    start1 = time.time()
    #try all the different rerootings and record the ones and their solutions with the best scores
    # print("The number of rerootings is %d" %len(all_rootings))
    for rooting in all_rootings:
        geneTree=rooting
        # start2 = time.time()
        MPR,cost=DTLOR_DP.DP(speciesTree, geneTree, phi, locusMapD, D, T, L, O, R)
        # end2 = time.time()
        # print("The time took to do one reconciliation is: %.4f" %(end2 - start2))
        if cost<best_score: 
            #if the score is better than current best
            #update best score, clear record and add new record
            best_score=cost
            bestMPRs=[]
            bestMPRs.append((geneTree,MPR))
        elif cost==best_score:
            bestMPRs.append((geneTree,MPR))
    #sample one MPR from the MPRs for this specific unrooted tree
    end1 = time.time()
    # print("The time took to reconcile all the rerootings is: %.4f" %(end1 - start1))
    optGeneRooting,optMPR=random.choice(bestMPRs) 

    return optGeneRooting,optMPR    

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


def readFamilies(familyFN,tree,genesO):
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
        mrca = L[1]

        geneTree = eval(L[2])
        recon = eval(L[3].replace('inf', "float('inf')")) # a hack!
        # eval didn't like the string inf.
                     
        lfL = L[4:]

        familiesO.initializeFamily(famNum,mrca,geneTree,recon)

        for lfStr in lfL:
            lfSplitL = lfStr.rstrip().split(',')
            locusFamNum=int(lfSplitL[0])
            lfMrca = lfSplitL[1]
            locusNum = eval(lfSplitL[2])
            geneL=[]
            for geneName in lfSplitL[3:]:
                geneNum = int(geneName.split('_')[0])
                geneL.append(geneNum)
            lfO = LocusFamily(famNum,locusFamNum,lfMrca,locusNum)
            lfO.addGenes(geneL,genesO)
            familiesO.addLocusFamily(lfO)

    f.close()
    return familiesO
