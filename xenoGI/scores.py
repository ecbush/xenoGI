import parasail,pickle,glob,statistics,numpy
from scipy.signal import find_peaks
from multiprocessing import set_start_method, Pool
from . import genomes
from . import trees
from . import Score

#### Functions


## raw similarity scores

def calcRawScores(fastaFilePath,numThreads,geneNames,gapOpen,gapExtend,matrix,scoresO):
    '''Get a global alignment based raw score for every edge in scoresO.'''

    # load sequences
    protFnL=glob.glob(fastaFilePath)
    seqD=genomes.loadProt(protFnL)
                
    # make list of sets of arguments to be passed to p.map. There
    # should be numThreads sets.
    argumentL = [([],seqD,gapOpen, gapExtend, matrix) for i in range(numThreads)]

    i=0
    for g1,g2 in scoresO.iterateEdgesByEndNodes():
        edgeNum = scoresO.endNodesToEdge(g1,g2)
        edgeT = edgeNum,geneNames.numToName(g1),geneNames.numToName(g2)
        argumentL[i%numThreads][0].append(edgeT)
        i+=1
        
    # run
    p=Pool(numThreads)
    scoresLL = p.map(rawScoreGroup, argumentL)
    p.close()
    p.join()

    
    # store in scoresO
    for scoresL in scoresLL:
        for edgeNum,sc in scoresL:
            scoresO.addScoreByEdge(edgeNum,sc,'rawSc')

    return scoresO

def rawScoreGroup(argT):
    '''Given a dictionary of sequences and a list of gene pairs, go
through each pair and get a needleman wunch based score.
    '''
    edgeL,seqD,gapOpen,gapExtend,matrix = argT
    matrix = eval(matrix) # turn it from string to parasail matrix object
    scoresL = []
    for edgeNum,g1,g2 in edgeL:
        scaled = rawScore(seqD[g1],seqD[g2],gapOpen,gapExtend,matrix)
        scoresL.append((edgeNum,scaled))
    return scoresL
      
   
def rawScore(s1,s2,gapOpen, gapExtend, matrix):
    '''Calculate score between a pair of protein sequences, based on a
global alignment. We scale the alignment score to be between 0 and 1,
based on the max and min possible scores for these sequences.'''

    # note. Here we've got a score 0-1 based only on the alignment
    # score, and not requiring the alignment itself. Right now
    # parasail doens't give the alignment. It might in the future, and
    # we could revisit this.
    
    r_s1s2 = parasail.nw_scan(s1,s2, gapOpen, gapExtend, matrix)

    if len(s1) < len(s2):
        r_self = parasail.nw_scan(s1,s1, gapOpen, gapExtend, matrix)
    else:
        r_self = parasail.nw_scan(s2,s2, gapOpen, gapExtend, matrix)

    sc = r_s1s2.score
    mx = r_self.score # max possible is shorter seq against itself.
    
    # lowest possible score, if we have gaps opposite all residues and
    # two opens. Note parasail does not count gap extend for the
    # residue where a gap is opened, hence the -2 in the extend
    # formula below.
    mn = - ( 2 * gapOpen + ( (len(s1)+len(s2) -2 ) * gapExtend ) )
    scaled = (sc - mn) / (mx - mn)
    
    return scaled



## normalized scores


def calcNormScores(tree,strainNum2StrD,blastFilePath,evalueThresh,scoresO,geneNames,aabrhFN):
    '''Given directory of blast output and a graph of raw similarity
scores, calculate normalized similarity scores by comparing each score
with the expection for two genes evolving at an "average" rate. We
obtain this expectation by finding the strains for our two genes, and
then getting all the scores from all around best reciprocal hits in
that pair of strains. For genes coming from different strains, we
normalize using the mean and std for aabrh core gene raw scores
between those strains. For genes from the same strain, we find the
nearest neigbhor strain(s) to that strain, and then get the mean and
std of scores between the starting strain and these. We then use this
to calculate the norm score. This approach to genes from the same
strain allows us to have norm scores for every pair of genes where
there is a raw score.

Here are some things to think about in considering this. The idea of a
norm score is that two genes who have evolved at an "average" rate
since their common ancestor should get a 0. In the case of genes from
different strains, it makes sense to normalize with known homologs in
the two. In genes of the same strain, the mrca is a duplication, and
could either be before or after the time of divergence of the nearest
neighbor. If it is before, this will tend to result in slightly
negative norm scores for genes that evolved at the "average" rate. If
the duplication comes after the most recent divergence, it will make
norm scores that are slightly positive for an average evolving pair of
genes. However, these biases are small. We are doing it this way
because it is an advantage (and a simplification) to have the same
sorts of scores available for all pairs of genes.
    '''

    strainNamesL=sorted([strainNum2StrD[leaf] for leaf in trees.leafList(tree)])
    aabrhL = createAabrhL(blastFilePath,strainNamesL,evalueThresh,aabrhFN)

    strainNumsL=sorted(leaf for leaf in trees.leafList(tree))
    scoresO.createAabrhScoreSummaryD(strainNumsL,aabrhL,geneNames)
   
    # loop over each edge in scoresO, normalizing score and saving there
    for gn1,gn2 in scoresO.iterateEdgesByEndNodes():
        sp1 = geneNames.numToStrainNum(gn1)
        sp2 = geneNames.numToStrainNum(gn2)
        rawSc=scoresO.getScoreByEndNodes(gn1,gn2,'rawSc')
            
        if sp1 == sp2:
            # find mean,std from scoresO.scoreSummaryD for the nearest
            # neighbors to this species.
            mean,std = getNearestNeighborAverageScore(sp1,tree,scoresO)
            normSc = normScore(rawSc,mean,std)
        else:
            # find mean,std from scoresO.scoreSummaryD.
            mean,std = scoresO.scoreSummaryD[(sp1,sp2)]
            normSc = normScore(rawSc,mean,std)
            
        scoresO.addScoreByEndNodes(gn1,gn2,normSc,'normSc')

    return scoresO,aabrhL


def createAabrhL(blastFilePath,strainNamesL,evalueThresh,aabrhFN):
    '''Get the sets of all around best reciprocal hits.'''

    blastDir = blastFilePath.split("*")[0]
    rHitsL=getAllReciprocalHits(blastDir,strainNamesL,evalueThresh)
    
    # get sets of genes in first species that have a reciprocal best
    # hit in each other species. this is to save time in next step.
    possibleGenesD = getPossibleGenesD(rHitsL)

    # then, out of the possible genes, get those that are all
    # pairwise reciprocal bets hits
    aabrhL = getAllAroundBRHL(possibleGenesD, rHitsL)

    # write it to file
    f=open(aabrhFN,'w')
    for orthoT in aabrhL:
        f.write("\t".join(orthoT)+"\n")
    f.close()
        
    return aabrhL


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

def loadOrthos(aabrhFN):
    '''Reads the all around best reciprocal hits orthologs file. One set
per line.'''
    f= open(aabrhFN,'r')

    orthoL=[]
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.rstrip().split('\t')
        orthoL.append(tuple(L))
            
    f.close()
    return orthoL

def getAllReciprocalHits(blastDir,strainNamesL,evalueThresh):
    '''return an upper-diagonal (N-1)xN matrix where each entry
    [i][j] (j > i) contains a dictionary of best reciprocal hits
    between species i and species j, as indexed in list strainNamesL;
    key is always gene in species i, and value gene in species j'''

    rHitsL = []
    
    for i in range(len(strainNamesL)): 
        rHitsL.append([])
        # fill in the rest of the row with the dictionary of best
        # reciprocal hits between species i and j (keyed by species i)
        for j in range(len(strainNamesL)):
            if j > i:
                rHitsL[i].append(getReciprocalHits(strainNamesL[i], strainNamesL[j], blastDir, evalueThresh))
            else:
                rHitsL[i].append(None)
    return rHitsL


def getReciprocalHits(strainName1, strainName2, blastDir, evalueThresh):
    '''Given strain names and blast file directory name, load hits between
two strains in each direction, then go through and keep only the
reciprocal hits. Returns dictionary where keys are genes in strain 1
and values are corresponding genes in strain 2.'''
    # get best hits of 1 blasted against 2...
    hits1D = getHits(blastDir+strainName1+'-'+strainName2+'.out', evalueThresh)
    # ...and 2 blasted against 1
    hits2D = getHits(blastDir+strainName2+'-'+strainName1+'.out', evalueThresh)

    # then store only the reciprocal best hits
    recipHitsD = {}
    for query in hits1D:
        hit = hits1D[query]
        nextHit = hits2D.get(hit)
        if nextHit == query:
            recipHitsD[query] = hit

    return recipHitsD

def getHits(fileName, evalueThresh):
    """Given a BLAST output file, returns a dictionary keyed by the genes
    in the query species, with the values being the top hit (if any)
    for those genes. Assumes the blast hits for each query are given
    from most to least significant. which appears to be the
    case. Thresholds for minimum similarity and maximum length
    difference are globally defined.
    """

    f = open(fileName, 'r')
    hitsD = {}

    queryGene = ""

    # read through entire file
    line = f.readline()

    while line != '':
        # break the line up into its components
        L = line.split('\t')

        if len(L)<12:
            # its not an info line (likely header)
            line = f.readline()
            continue

        # gene names come first
        queryGene = L[0]
        hit = L[1]
        evalue = float(L[10])

        if evalue < evalueThresh:
            # store the query and hit if they meet thresholds
            hitsD[queryGene] = hit

        # we only want the first hit for any query gene. Change this
        # later to check to get best score even if doesn't come first.
        while line.split('\t')[0] == queryGene:
            line = f.readline()

    f.close()
    return hitsD

def getPossibleGenesD(rHitsL):
    '''return dict containing an entry for each gene in species 0 that has
    a best reciprocal hit in all other species.  format of an entry is
    gene:(tuple of recip hit genes in other species). We do this to
    quickly eliminate things we aren't interested in.
    '''
    possibleGenesD = {}
    firstRow = rHitsL[0]
    firstDict = firstRow[1]
    # we just look at first dict, since the all around brh sets will
    # have a gene in this with matches in all strains.
    for gene in firstDict:
        isPossible = True
        hitsL = [firstDict[gene]]
        for colInd in range(2, len(firstRow)):
            compareDict = firstRow[colInd]
            hit = compareDict.get(gene)
            if hit == None:
                isPossible = False
                break
            else:
                hitsL.append(hit)
        if isPossible:
            possibleGenesD[gene] = tuple(hitsL)
    return possibleGenesD


def getAllAroundBRHL(possibleGenesD, rHitsL):
    '''Given a matrix of reciprocal hit dictionaries in rHitsL, and a dict
of possible hits, return a list of all around best reciprocal hits. In
the form of a list of lists, where each sublist is one set of
orthologs.'''
    
    aabrhL = []

    genesL = list(possibleGenesD.keys())
    genesL.sort() # sort so we get consistent order
    
    for gene in genesL:
        hitsT = possibleGenesD[gene]

        # check that each gene in the list of hitsT is best reciprocal 
        # hit of all following genes
        isOrth = True
        for hit1Index in range(len(hitsT)):
            hit1 = hitsT[hit1Index]

            # check that current gene in list is best reciprocal hit
            # of all following genes in list
            for hit2Index in range(hit1Index+1, len(hitsT)): # j always > i
                hit2 = hitsT[hit2Index]
                sp1sp2Dict = rHitsL[hit1Index+1][hit2Index+1]
                match = sp1sp2Dict.get(hit1)
                if match != hit2: # oh noes! not reciprocal!
                    isOrth = False
                    break

            # if any pair not reciprocal best hits, then we fail
            if not isOrth:
                break

        # if we get through all genes without finding that some aren't
        # reciprocal best hits, then we keep this list of genes
        if isOrth:
            aabrhL.append((gene,)+hitsT)
    return aabrhL

def normScore(rawScore,mean,std):
    '''Given a raw score, mean and std, return a
score normalized by std and centered around zero.'''
    norm = ( rawScore - mean ) / std
    return norm

## synteny scores

def calcSynScores(scoresO,geneNames,geneOrderT,synWSize,tree,numSynToTake,numThreads):
    '''Calculate the synteny score between two genes and add to edge
attributes of scoresO. We only bother making synteny scores for those
genes that have an edge in scoresO.
    '''
    
    neighborTL = createNeighborL(geneNames,geneOrderT,synWSize)

    # make list of groups of arguments to be passed to p.map. There
    # should be numThreads groups.
    argumentL = [([],tree,neighborTL,numSynToTake,geneNames,scoresO) for i in range(numThreads)]

    i=0
    for gn1,gn2 in scoresO.iterateEdgesByEndNodes():
        argumentL[i%numThreads][0].append((gn1,gn2))
        i+=1

    p=Pool(numThreads) # num threads
    synScoresLL = p.map(synScoreGroup, argumentL)
    p.close()
    p.join()
    
    # add to scores object
    for synScoresL in synScoresLL:
        for gn1,gn2,sc in synScoresL:
            scoresO.addScoreByEndNodes(gn1,gn2,sc,'synSc')

    return scoresO

def createNeighborL(geneNames,geneOrderT,synWSize):
    '''Return a list which specifies the neighbors of each gene. Index of
list corresponds to gene number, and the value located at that index
is a tuple of all genes within a synWSize window. e.g. synWSize 10 means we
go 5 genes in either direction.'''

    lenInEitherDirec = int(synWSize/2)
    
    neighborTL = [None for x in geneNames.nums]
    
    for contigT in geneOrderT:
        if not contigT == None:
            for geneNumT in contigT:
                for i in range(len(geneNumT)):
                    end = i + lenInEitherDirec
                    st = i-lenInEitherDirec if i-lenInEitherDirec>0 else 0 # st can't be less than 0
                    L = list(geneNumT[st:end])
                    L.remove(geneNumT[i])
                    neighborTL[geneNumT[i]] = tuple(L)

    return neighborTL

def synScoreGroup(argsT):
    '''Given an argument list, including pairs of genes, calculate synteny
    scores. This function is intended to be called by p.map.
    '''

    edgeL,tree,neighborTL,numSynToTake,geneNames,scoresO = argsT

    outL=[]
    for gn1,gn2 in edgeL:
        outL.append(synScore(scoresO,gn1,gn2,tree,neighborTL,numSynToTake,geneNames))
        
    return outL
        
def synScore(scoresO,gn1,gn2,tree,neighborTL,numSynToTake,geneNames):
    '''Given two genes, calculate a synteny score for them. We are given
    the genes, neighborTL, which contains lists of neighbors for each
    gene. For the two sets of neighbors, we find the numSynToTake top
    pairs, and return the average of their scores. The approach is
    greedy. We find the pair with the best score, add it, then remove
    those genes and iterate.
    '''

    # get the min possible score for these two species (this is
    # really for the case of using normalized scores, where it
    # varies by species pair.)
    sp1 = geneNames.numToStrainNum(gn1)
    sp2 = geneNames.numToStrainNum(gn2)

    if sp1 == sp2:
        # for same species, norm scores based on comparison with
        # nearest neighbor
        mean,std = getNearestNeighborAverageScore(sp1,tree,scoresO)
    else:
        mean,std = scoresO.scoreSummaryD[(sp1,sp2)]

    minNormScore = normScore(0,mean,std) # min raw score is 0
    L1 = list(neighborTL[gn1])
    L2 = list(neighborTL[gn2])

    topScL= [minNormScore] * numSynToTake

    for i in range(numSynToTake):
        ind1,ind2,sc = topScore(L1,L2,scoresO)
        if sc == -float('inf'):
            break
        topScL[i] = sc
        del(L1[ind1]) # remove this index
        del(L2[ind2])
        
    synSc = sum(topScL) / numSynToTake
    
    return gn1, gn2, synSc

def topScore(L1,L2,scoresO):
    '''Find the best norm score between genes in L1 and L2. Return the index of
each and the score.'''
    besti1 = 0
    besti2 = 0
    bestSc = -float('inf')

    for i1,gn1 in enumerate(L1):
        for i2,gn2 in enumerate(L2):
            if scoresO.isEdgePresentByEndNodes(gn1,gn2):
                sc = scoresO.getScoreByEndNodes(gn1,gn2,'normSc')
                if sc > bestSc:
                    bestSc = sc
                    besti1 = i1
                    besti2 = i2
    return besti1,besti2,bestSc


## Core synteny scores

def calcCoreSynScores(scoresO,aabrhL,geneNames,geneOrderT,coreSynWsize):
    '''Calculate synteny scores based on core genes given in
aabrhL. Scores are between 0 and 1, giving the percentage of syntenic
genes shared.'''

    geneToAabrhT = createGeneToAabrhT(aabrhL,geneNames)
    coreSyntenyT = createCoreSyntenyT(geneToAabrhT,geneOrderT,coreSynWsize)

    # Not parallelized. Pretty fast already.

    # loop over all edges in scoresO, adding corresponding core syn score
    for gn1,gn2 in scoresO.iterateEdgesByEndNodes():
        coreSynSc=coreSynScore(coreSyntenyT[gn1],coreSyntenyT[gn2],coreSynWsize)
        scoresO.addScoreByEndNodes(gn1,gn2,coreSynSc,'coreSynSc')
        
    return scoresO

def createGeneToAabrhT(aabrhL,geneNames):
    '''Create a tuple where the index corresponds to gene number and the
value at that location is the number of the aabrh group to which the
gene belongs, or None. The aabrh groups are numbered, simply based on
the index where they occur in aabrhL.'''

    geneToAabrhL =[None] * len(geneNames.nums)

    for aabrhNum in range(len(aabrhL)):
        for geneName in aabrhL[aabrhNum]:
            geneNum = geneNames.nameToNum(geneName)
            geneToAabrhL[geneNum] = aabrhNum
    return tuple(geneToAabrhL)

def createCoreSyntenyT(geneToAabrhT,geneOrderT,coreSynWsize):
    '''Create and return core synteny tuple. The index of this corresponds
to gene number. The value at that index is a tuple of the aabrh
numbers for the syntenic genes within coreSynWsize.'''

    coreSyntenyL =[None] * len(geneToAabrhT)
    geneAabrhOrderL = createGeneAabrhOrderL(geneToAabrhT,geneOrderT)
    for element in geneAabrhOrderL:
        if element != None:
            # internal nodes are None, having no genes themselves
            for geneNumT,aabrhNumT in element:
                for pos in range(len(geneNumT)):
                    geneNum=geneNumT[pos]
                    coreSyntenyL[geneNum]=getAabrhContext(pos,aabrhNumT,coreSynWsize)
    return tuple(coreSyntenyL)
                    
def getAabrhContext(pos,aabrhNumT,coreSynWsize):
    '''Starting at index pos in aabrhNumT, return the coreSynWsize/2 aabrh
group numbers to each side (not including the gene at pos).'''

    numToTakeEachSide = int(coreSynWsize/2)
    coreSynL=[]

    # go forward from gene
    ct=0
    j=pos+1
    while j < len(aabrhNumT) and ct < numToTakeEachSide:
        if aabrhNumT[j]!=None:
            coreSynL.append(aabrhNumT[j])
            ct+=1
        j+=1

    # go backward from gene
    ct=0
    j=pos-1
    while j >= 0 and ct < numToTakeEachSide:
        if aabrhNumT[j]!=None:
            coreSynL.append(aabrhNumT[j])
            ct+=1
        j-=1

    return tuple(coreSynL)
            
def createGeneAabrhOrderL(geneToAabrhT,geneOrderT):
    '''Make geneAabrhOrderL, which follows the structure of geneOrderT,
but instead of only having geneNumT's, it has (geneNumT,aabrhNumT)
pairs at its lowest level.
    '''
    geneAabrhOrderL=[]
    for contigT in geneOrderT:
        if contigT == None:
            geneAabrhOrderL.append(None)
        else:
            contigL=[]
            for geneNumT in contigT:
                aabrhNumT = convertGeneToAabrhTuple(geneNumT,geneToAabrhT)
                contigL.append((geneNumT,aabrhNumT))
            geneAabrhOrderL.append(tuple(contigL))
    return geneAabrhOrderL

def convertGeneToAabrhTuple(geneNumT,geneToAabrhT):
    '''Take a tuple of gene numbers and make a corresponding tuple where
every gene number is replaced by its corresponding aabrh number, or
None if it has none.'''
    aabrhNumL=[]
    for geneNum in geneNumT:
        aabrhNumL.append(geneToAabrhT[geneNum])
    return tuple(aabrhNumL)

def coreSynScore(synT1,synT2,coreSynWsize):
    '''Given two tuples, determine the number of shared core genes in
each. Return this, divided by the core synteny window size, which is
the number of genes we consider around each gene. Values will be
between 0 and 1.'''
    ct=0
    for cgene in synT1:
        if cgene in synT2:
            ct+=1
    return ct/coreSynWsize

## Histograms and thresholds

def getHomologyRawThresholdD(scoresO):
    '''For each pair of strains (including vs. self) determine a minimum
raw score above which we take scores to indicate homology. Return in a
dictionary keyed by strain pair. It also returns a boolean flag,
homologousPeakMissing. This will be true if any of the strain pairs
failed to yield a peak on the right side of the histogram
corresponding to homologous scores. This function works as
follows. For a given strain pair, it expects the histogram of scores
to have two peaks. A left non homologous peak, and a right homologous
peak. It attempts to find the right peak, and then finds the left
peak. It takes the right base of the left non-homologous peak to be
the threshold. Unless the left base of the right peak is actually smaller, in which case it takes that.

    '''

    defaultHomologyThreshold = 0.75
    
    numBins = 80
    binWidth = 1.0/numBins # since scores range from 0-1

    homologousPeakMissing = False
    homologyRawThresholdD = {}
    for strainPair in scoresO.getStrainPairs():
        scoreIterator = scoresO.iterateScoreByStrainPair(strainPair,'rawSc')
        binHeightL,indexToBinCenterL = scoreHist(scoreIterator,numBins)
        homologPeakLeftExtremePos=homologPeakChecker(binHeightL,indexToBinCenterL,binWidth)

        if homologPeakLeftExtremePos == float('inf'):
            homologousPeakMissing = True
            # if we can't find the homologous peak, there's really no
            # point in examining the non-homologous peak. Just use
            # default threshold.
            threshold = defaultHomologyThreshold
        else:
            threshold = getMinRawThreshold(binHeightL,indexToBinCenterL,binWidth,homologPeakLeftExtremePos)
        homologyRawThresholdD[strainPair] = threshold


    return homologousPeakMissing,homologyRawThresholdD


def scoreHist(scoreIterator,numBins):
    '''Get a histogram with numpy, and return the bin height, and also a
list of indices to the middle position of each bin (in terms of the x value).'''
    binHeightL,edges = numpy.histogram(list(scoreIterator),bins=numBins,density=True)

    # make a list where the indices correspond to those of binHeightL,
    # and the values give the score value at the center of that bin
    indexToBinCenterL = []
    for i in range(1,len(edges)):
        center = (edges[i]+edges[i-1])/2
        indexToBinCenterL.append(center)

    return binHeightL,indexToBinCenterL

def homologPeakChecker(binHeightL,indexToBinCenterL,binWidth):
    '''Function to check for a peak due to homogology (right peak in
histogram). If such a peak exists, this function returns the position
(in score units) of the left most base of that peak. If no such peak
exits, this function returns infinity.

    '''
    peakL = [] # collect them all first
    
    # case 1 (normal case)
    homologPeakWidth = 0.10
    widthRelHeight = 0.9
    homologRequiredProminence = 0.2
    homologLeftPeakLimit = 0.65
    homologRightPeakLimit = 1.0
    
    L = findPeaksOneCase(binHeightL,indexToBinCenterL,binWidth,homologPeakWidth,widthRelHeight,homologRequiredProminence,homologLeftPeakLimit,homologRightPeakLimit)
    peakL.extend(L)

    # case 2 (extreme prominence. But allow to be very narrow)
    homologPeakWidth = 0
    widthRelHeight = 0.9
    homologRequiredProminence = 6
    homologLeftPeakLimit = 0.90
    homologRightPeakLimit = 1.0
    
    # in order to get a rightmost peak, if any, we add a dummy bin of height 0 on right
    # add to indexToBinCenterL for case of right base of a peak in the last bin.
    tempBinHeightL = numpy.append(binHeightL,0)
    tempIndexToBinCenterL = numpy.append(indexToBinCenterL,1)
    L = findPeaksOneCase(tempBinHeightL,tempIndexToBinCenterL,binWidth,homologPeakWidth,widthRelHeight,homologRequiredProminence,homologLeftPeakLimit,homologRightPeakLimit)
    peakL.extend(L)

    # case 3 (wide width with low prominence)
    homologPeakWidth = 0.25
    widthRelHeight = 0.9
    homologRequiredProminence = 0.10
    homologLeftPeakLimit = 0.65
    homologRightPeakLimit = 1.0
    
    L = findPeaksOneCase(binHeightL,indexToBinCenterL,binWidth,homologPeakWidth,widthRelHeight,homologRequiredProminence,homologLeftPeakLimit,homologRightPeakLimit)
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

def getMinRawThreshold(binHeightL,indexToBinCenterL,binWidth,homologPeakLeftExtremePos):
    '''Given a list of bin heights and another list giving the score
values at the middle of each bin, determine a threshold above which a
score should be taken to indicate homology. This function assumes
there is a right homologous peak present and takes the left extreme of
this peak as input.

    '''

    nonHomologPeakWidth = 0.15
    widthRelHeight = 0.9
    nonHomologPeakProminence = 2
    defaultLeftBasePos = 0.8
    nonHomologLeftPeakLimit = 0
    nonHomologRightPeakLimit = 0.6

    L = findPeaksOneCase(binHeightL,indexToBinCenterL,binWidth,nonHomologPeakWidth,widthRelHeight,nonHomologPeakProminence,nonHomologLeftPeakLimit,nonHomologRightPeakLimit)

    L.sort(reverse=True) # in the unlikely case there's more than one
    peakHeight,peakPos,leftExtremeOfPeakPos,rightExtremeOfPeakPos = L[0]

    # we now find the minimum of these two. We're after a threshold
    # where we're confident that things below it are definitely not
    # homologous.
    threshold = min(rightExtremeOfPeakPos,homologPeakLeftExtremePos)

    return threshold
    
## Graph I/O

def writeScores(scoresO,geneNames,scoresFN):
    '''Write graph G to file. If scoresFN has the .bout extension, write
binary pickle of G, otherwise write in text output format.'''

    if scoresFN.split('.')[-1] == 'bout':
        scoresO.writeScoresBinary(scoresFN)
    else:
        scoresO.writeScoresText(geneNames,scoresFN)


def readScores(scoresFN,geneNames=None):
    '''Read scores from file creating a Score object of scores. If
scoresFN has the .bout extension, read binary pickle of object,
otherwise read text format.'''
    if scoresFN.split('.')[-1] == 'bout':
        scoresO = Score.Score.readScoresBinary(scoresFN)
    else:
        scoresO = Score.Score.readScoresText(scoresFN,geneNames)
    return scoresO
