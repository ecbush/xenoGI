import sys,numpy,os,random,glob,copy
from scipy.signal import find_peaks
from Bio import Phylo
from multiprocessing import Pool
from . import trees,scores,DTLOR_DP
from .Family import *
from .Tree import *
from .analysis import printTable
#import matplotlib.pyplot as plt

#### Main function

def createFamiliesO(speciesRtreeO,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,outputSummaryF):
    '''Main function to create families. First creates an initial families
object, then reconciles those against the species tree to make an
originFamilies object.
    '''

    # threshold for maximum size of gene tree, if above, split
    maxFamilySize=25

    # define some variables
    initFamilyFN = paramD['initFamilyFN']
    originFamilyFN =  paramD['originFamilyFN']
    geneInfoFN = paramD['geneInfoFN']
    iFamGeneTreeFileStem = paramD['iFamGeneTreeFileStem']

    # checks
    homologyCheck(genesO,aabrhHardCoreL,scoresO,outputSummaryF,paramD)

    # create initial families
    initialFamiliesO, locusMapD=createInitialFamiliesO(maxFamilySize,speciesRtreeO,scoresO,genesO,aabrhHardCoreL,paramD)

    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT)) # needed for writeFamilies
    writeFamilies(initialFamiliesO,initFamilyFN,genesO,strainNamesT,paramD)
    print("Initial families:",file=outputSummaryF)
    writeFamilyFormationSummary(initialFamiliesO,outputSummaryF)
    print("Initial families written out to %s"%initFamilyFN,file=sys.stderr)

    # make gene family trees
    # TEMP COMMENT
    trees.makeGeneFamilyTrees(paramD,genesO,initialFamiliesO,iFamGeneTreeFileStem) #create gene tree for each initial family 
    print("Finished making gene trees")

    # load gene trees
    singleGeneInitFamNumL,multifurcatingL,bifurcatingL = loadGeneTrees(paramD,initialFamiliesO,iFamGeneTreeFileStem)

    # strip out trees from multifurcatingL (until reconciliation can handle them)
    multifurcatingL = [initFamNum for initFamNum,_ in multifurcatingL]
    
    # reconcile

    # TEMP
    with open("reconTemp.tsv","w") as f:
        print("",file=f) # blank temp file
    # TEMP
    
    initialFamiliesO = reconcileAllGeneTrees(speciesRtreeO,bifurcatingL,initialFamiliesO,locusMapD,genesO,paramD)
    writeFamilies(initialFamiliesO,initFamilyFN,genesO,strainNamesT,paramD)
    print("Initial families updated with reconciliations and gene trees",file=sys.stderr)    

    # create origin families
    originFamiliesO = createOriginFamiliesO(speciesRtreeO,singleGeneInitFamNumL,multifurcatingL,bifurcatingL,initialFamiliesO,genesO)

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

def createInitialFamiliesO(maxFamilySize,speciesRtreeO,scoresO,genesO,aabrhHardCoreL,paramD):
    '''Given a scoresO object, create initial genes families based on
    blast connectivity and return as a familiesO object.
    '''
    
    # get synteny thesholds for locus family formation
    synThresholdD = getSynThresholdD(paramD,scoresO,genesO,aabrhHardCoreL,speciesRtreeO)

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
    initFamiliesO = Families(speciesRtreeO)

    # store in familiesO object
    # these locusFam Ids need to be unique across all, start counting from 1. EB: WHY not 0?
    famNumCounter=1
    locusFamNumCounter=1
    totalAddedTolocusFamilies=0
    for size, familyS, locusFamLL in initFamilyDivIntoLocFamsL:
        locusFamNumCounter, famNumCounter,totalAddedTolocusFamilies=addFamilyToInitialFamiliesO(familyS,initFamiliesO,famNumCounter,locusFamLL,locusMapD,locusFamNumCounter,genesO,speciesRtreeO,totalAddedTolocusFamilies)
    return initFamiliesO, locusMapD

def getSynThresholdD(paramD,scoresO,genesO,aabrhHardCoreL,speciesRtreeO):
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

def addFamilyToInitialFamiliesO(familyS,initFamiliesO,famNumCounter,locusFamLL,locusMapD,locusFamNumCounter,genesO,speciesRtreeO,totalAddedTolocusFamilies):
    '''Given a families object, a family set and a list of locus families
in that family, add the family to the object.
    '''
    speciesL=list(set([genesO.numToStrainName(gene) for gene in familyS]))
    mrca=speciesRtreeO.findMrca(speciesL)
    #add each initial family as a Family object (still empty)
    initFamiliesO.initializeFamily(famNumCounter,mrca) 
    for locusFamilyL in locusFamLL: #locusFamilyL contains all the genes in that lf
        speciesL=[]
        for gene in locusFamilyL:
            locusMapD[gene]=locusFamNumCounter
            speciesL.append(genesO.numToStrainName(gene))
        lfMrca=speciesRtreeO.findMrca(speciesL)
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

def loadGeneTrees(paramD,initialFamiliesO,iFamGeneTreeFileStem):
    '''Load gene trees, returning 3 things: A list of initial family
numbers from single gene families, a list of initial family numbers
from families with multifurcations, and a list of trees from purely
binary families with more than 1 genes.
    '''

    geneFamilyTreesDir = paramD['geneFamilyTreesDir']

    # dir of gene trees lacks single gene families. Get their numbers.
    singleGeneInitFamNumL=[]
    for initFam in initialFamiliesO.iterFamilies():
        if len(tuple(initFam.iterGenes())) == 1:
            singleGeneInitFamNumL.append(initFam.famNum)

    # load gene trees, divide into bifurcating vs. multifurcating
    allGtFilePath = os.path.join(geneFamilyTreesDir,iFamGeneTreeFileStem+'*.tre')
    allTreeFN_L=list(sorted(glob.glob(allGtFilePath)))

    multifurcatingL = []
    bifurcatingL = []
    for geneTreeFN in allTreeFN_L:
        geneUtreeO = Utree()
        geneUtreeO.fromNewickFile(geneTreeFN)
        initFamNum = int(geneTreeFN.split(iFamGeneTreeFileStem)[1].split('.tre')[0].lstrip('0'))
        if geneUtreeO.multifurcatingNodes() != []:
            multifurcatingL.append((initFamNum,geneUtreeO))
        else:
            bifurcatingL.append((initFamNum,geneUtreeO))
        
    return singleGeneInitFamNumL,multifurcatingL,bifurcatingL

def reconcileAllGeneTrees(speciesRtreeO,geneTreeL,initialFamiliesO,locusMapD,genesO,paramD):
    '''Reconcile gene family trees to the species tree using the DTLOR algorithm.'''

    D=float(paramD["duplicationCost"])
    T=float(paramD["transferCost"])
    L=float(paramD["lossCost"])
    O=float(paramD["originCost"])
    R=float(paramD["rearrangeCost"])

    argumentL = []
    for initFamNum,geneUtreeO in geneTreeL:

        # in loop
        tipMapD=getTipMapping(geneUtreeO,genesO)

        # Make new gtLocusMapD with only those genes in geneUtreeO
        gtLocusMapD = reduceLocusMap(geneUtreeO,locusMapD)
        
        # add to argumentL
        argT = (initFamNum,speciesRtreeO,geneUtreeO,tipMapD,gtLocusMapD,D,T,L,O,R)
        argumentL.append(argT)

    # run on multiple processors
    with Pool(processes=paramD['numProcesses']) as p:
        for initFamNum,optGeneRtreeO,optMPR,minCost in p.imap_unordered(reconcile, argumentL):
            
            reconD = convertReconBranchToNode(optMPR)
            
            # sanity check
            #costCheck(minCost,reconD,D,T,L,O,R)

            # TEMP
            sm=costSum(reconD,D,T,L,O,R)
            if round(sm,3) != round(minCost,3):
                print("Cost problem, initfam:",initFamNum,"cost",minCost)

            # TEMP
            def findArg(initFamNum,argumentL):
                for argT in argumentL:
                    if argT[0] == initFamNum:
                        return argT
                return None

            argT=findArg(initFamNum,argumentL)
            with open("reconTemp.tsv","a") as f:
                outL = [str(initFamNum),optGeneRtreeO.fileStr(),str(optMPR),str(minCost),str(argT)]
                outStr = "\t".join(outL)
                print(outStr,file=f)
                
            # TEMP
                
            # store
            ifam = initialFamiliesO.getFamily(initFamNum)
            ifam.addGeneTree(optGeneRtreeO)
            ifam.addReconciliation(reconD)

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
    bestMPRs=[]
    for geneRtreeO in geneUtreeO.iterAllRootedTrees():

        try:
            geneTreeD = geneRtreeO.createDtlorD(False) # put in dp format
        except:
            print(initFamNum)
            sys.exit()
            
        MPR,cost=DTLOR_DP.DP(speciesTreeD,geneTreeD,tipMapD,gtLocusMapD,D,T,L,O,R)
                
        if cost<minCost: 
            #if the score is better than current best
            #update best score, clear record and add new record
            minCost=cost
            bestMPRs=[]
            bestMPRs.append((geneRtreeO,MPR,minCost))
        elif cost==minCost:
            bestMPRs.append((geneRtreeO,MPR,minCost))

    #sample one MPR from the MPRs for this specific unrooted tree
    optGeneRtreeO,optMPR,minCost=random.choice(bestMPRs) 
        
    return initFamNum,optGeneRtreeO,optMPR,minCost

def costCheck(minCost,reconD,D,T,L,O,R):
    '''Sanity check to ensure the events in reconD sum to minCost.'''
    sm = costSum(reconD,D,T,L,O,R)
    assert round(sm,3) == round(minCost,3), "Events in this reconcilation don't sum to minCost."

def costSum(reconD,D,T,L,O,R):
    '''Sum the costs implied by the recon in reconD.'''
    sm = 0
    for valL in reconD.values():
        for eventType,_,_,_ in valL:
            if eventType == 'D':
                sm+=D
            elif eventType == 'T':
                sm+=T
            elif eventType == 'L':
                sm+=L
            elif eventType == 'O':
                sm+=O
            elif eventType == 'R':
                sm+=R
            # S, C and N events have no cost
    return sm
    
def createOriginFamiliesO(speciesRtreeO,singleGeneInitFamNumL,multifurcatingL,bifurcatingL,initialFamiliesO,genesO):
    '''Create and return an originFamilies object, based on the initial families and recocniliations.'''

    originFamiliesO = Families(speciesRtreeO)

    # add in families where we have no reconciliation
    originFamiliesO = addOriginFamilyFromInitialFamiliesO(singleGeneInitFamNumL,False,initialFamiliesO,originFamiliesO,genesO)
    originFamiliesO = addOriginFamilyFromInitialFamiliesO(multifurcatingL,True,initialFamiliesO,originFamiliesO,genesO)

    # add in families from reconciliation
    for initFamNum,geneTree in bifurcatingL:
        iFam = initialFamiliesO.getFamily(initFamNum)
        sourceFam = iFam.famNum
        originFamiliesO = addOriginFamilyFromReconciliation(iFam.geneRtreeO,iFam.reconD,originFamiliesO,sourceFam,genesO)

    return originFamiliesO
        
def addOriginFamilyFromInitialFamiliesO(initFamNumL,isMultiFurc,initialFamiliesO,originFamiliesO,genesO):
    '''For cases where there is no reconciliation, base origin family on
corresponding initial family.
    '''
    for initFamNum in initFamNumL:
        initFamO = initialFamiliesO.getFamily(initFamNum)
        famNum=originFamiliesO.getNumFamilies() # num for new family
        if isMultiFurc:
            originFamiliesO.initializeFamily(famNum,initFamO.mrca)
        else:
            # should be single gene family
            geneRtreeO,reconD = createSingleGeneFamilyGeneTreeRecon(initFamO,genesO)
            originFamiliesO.initializeFamily(famNum,initFamO.mrca,geneRtreeO,reconD)
            
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
    '''Given a rooted gene tree and a reconciliation (raw output from
DTLOR_DP, add origin families. One origin family for each origin
event.
    '''
    
    # get origin defined trees
    branchOriginL = getBranchesWithSpecifiedEvents(reconD,"O")
    originTreeL = splitTreeByOrigin(geneRtreeO,branchOriginL)
    
    # extract parts of reconD corresponding to each tree
    for geneRtreeO,splitReconD in getGeneTreeReconPairs(originTreeL,reconD):
        famNum = originFamiliesO.getNumFamilies()
        keyNB = (geneRtreeO.rootNode,'n')
        speciesMrca = splitReconD[keyNB][0][1] # remember, the value in reconD is a list of events
        originFamiliesO.initializeFamily(famNum,speciesMrca,geneRtreeO,splitReconD,sourceFam)
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
    
#### Convert reconciliation format

def convertReconBranchToNode(brReconD):
    '''Take a reconciliation dictionary in the branch based format. (The
format produced by the dtlor dp function, which represents the
placement of gene tree branches on the species tree). Return a
dictionary where reconciliation events are explicitly characterized as
falling at nodes or branches in the gene tree.

    '''

    rootKey = ''
    for key in brReconD:
        if key[0] == 'root':
            rootKey = key
            
    # get list of events. initial synLocAtTop should be '*'
    nodeBranchEventsL = getNodeBranchEventsL(brReconD,rootKey,'*')

    # create dict with keys (geneTreeLoc,'n/b')
    reconD = {}
    for eventKey,eventValue in nodeBranchEventsL:
        # eventKey is of form (geneTreeLoc,'n/b'). May not be unique
        
        if eventKey in reconD:
            reconD[eventKey].append(eventValue)
        else:
            reconD[eventKey] = [eventValue]
        
    return reconD

def getNodeBranchEventsL(brReconD,brKey,synLocAtTop):
    '''Traverse the branch based reconciliation dict in preorder starting
at brKey. synLocAtTop is the syntenic locus at the top (beginning) of
the branch brKey is on.

    '''
    if brKey not in brReconD:
        return [] # reached tip of gene tree
    else:
        _, child1Mapping, child2Mapping = brReconD[brKey]
        nodeBranchEventsL = parseOneBranchEntry(brReconD,brKey,synLocAtTop)
        synLocAtBottom = brKey[2]
        lL = getNodeBranchEventsL(brReconD,child1Mapping,synLocAtBottom) # left
        rL = getNodeBranchEventsL(brReconD,child2Mapping,synLocAtBottom) # right
        return nodeBranchEventsL + lL + rL

def parseOneBranchEntry(brReconD,brKey,synLocAtTop):
    '''Given a key to the branch based reconciliation dictionary, extract
the event or events (can have multiple when O or R occurs), and return
as a list of tuples. This list has all the info from this entry
(though there may be other entries which are also on the same gene
tree branch). Here is the form of the tuples

((geneTreeLoc,geneTreeNB),(eventType,speciesTreeLoc,spTreeNB,synLocAtBottom))

Each tuple corresponds to one event.

synLocAtTop is the syntenic locus at the top (beginning) of the branch
this entry is on. We use it to identify R events.
    '''
    nodeBranchEventsL = [] # for output
    
    geneTreeLoc,speciesTreeLoc,synLocAtBottom = brKey
    eventType, child1Mapping,child2Mapping = brReconD[brKey]

    # recognize N events as outside species tree
    if eventType == 'N':
        # we are outside the species tree
        speciesTreeLoc = 'xeno'
        
    # we want to refer to co-termination events with 'M', not 'C'
    if eventType == 'C':
        eventType = 'M'

    # determine branch or node for geneTreeLoc,speciesTreeLoc ('n' or 'b')
    geneTreeNB,spTreeNB = determineNodeOrBranch(eventType)

    eventKey = (geneTreeLoc,geneTreeNB)
    eventValue = (eventType,speciesTreeLoc,spTreeNB,synLocAtBottom)
    nodeBranchEventsL.append((eventKey,eventValue))

    # check for O,R events if eventType is D,T,S or M (not L, because
    # this may cause it to show up multiple times)
    if eventType in 'DTSM':
        if synLocAtTop=='*' and synLocAtBottom!='*':
            # O event
            eventKey = (geneTreeLoc,'b') # the branch leading to this DTSM event
            eventValue = ('O',speciesTreeLoc,'b',synLocAtBottom)
            nodeBranchEventsL.append((eventKey,eventValue))
        elif synLocAtTop!='*' and synLocAtTop != synLocAtBottom:
            # R event
            eventKey = (geneTreeLoc,'b') # the branch leading to this DTSM event
            eventValue = ('R',speciesTreeLoc,'b',synLocAtBottom)
            nodeBranchEventsL.append((eventKey,eventValue))
    
    return nodeBranchEventsL

def determineNodeOrBranch(eventType):
    '''Determine whether the gene and species locations are on nodes or
branches based on the event type.'''

    if eventType == 'D': # duplication
        # implies geneTreeLoc is a node, speciesTreeLoc is a branch
        return ('n','b')
    elif eventType == 'T': # transfer
        # implies geneTreeLoc is a node, speciesTreeLoc is a branch
        return ('n','b')
    elif eventType == 'L': # loss
        # implies geneTreeLoc is a branch, speciesTreeLoc is a node
        return ('b','n')
    elif eventType == 'O': # origin
        # implies geneTreeLoc is a branch, speciesTreeLoc is a branch
        return ('b','b')
    elif eventType == 'R': # rearrangement
        # implies geneTreeLoc is a branch, speciesTreeLoc is a branch
        return ('b','b')
    elif eventType == 'M': # cotermination
        # implies geneTreeLoc is a node, speciesTreeLoc is a node
        return ('n','n')
    elif eventType == 'S': # cospeciation
        # implies geneTreeLoc is a node, speciesTreeLoc is a node
        return ('n','n')
    if eventType == 'N': # null, when outside species tree
        # implies geneTreeLoc is a node, speciesTreeLoc is a branch
        return ('n','b')
    else:
        return # should never get here

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

def readFamilies(familyFN,speciesRtreeO,genesO):
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

        geneRtreeO = Rtree()
        geneRtreeO.fromString(L[2])
        
        recon = eval(L[3].replace('inf', "float('inf')")) # a hack!
        # eval didn't like the string inf.
        sourceFam = eval(L[4])
        familiesO.initializeFamily(famNum,mrca,geneRtreeO,recon,sourceFam)
        
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
