import parasail,glob,statistics,sys
from multiprocessing import set_start_method, Pool
from . import genomes
from . import trees
from . import Score


#### Global variables for use with multithreading

sharedScoresO = Score.sharedScore()

#### raw similarity scores

def calcRawScores(paramD,geneNames,scoresO):
    '''Get a global alignment based raw score for every edge in scoresO.'''

    numThreads = paramD['numThreads']
    gapOpen = paramD['gapOpen']
    gapExtend = paramD['gapExtend']
    matrix = paramD['matrix']
    
    # load sequences
    protFnL=glob.glob(paramD['fastaFilePath'])
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

    # run in multiple processes
    with Pool(processes=numThreads) as p:
        # store the results to scoresO as they come in
        for scoresL in p.imap_unordered(rawScoreGroup, argumentL):
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


#### synteny scores

def calcSynScores(scoresO,geneNames,geneOrderT,paramD,tree):
    '''Calculate the synteny score between two genes and add to edge
attributes of scoresO. We only bother making synteny scores for those
genes that have an edge in scoresO.
    '''

    synWSize = paramD['synWSize']
    numSynToTake = paramD['numSynToTake']
    numThreads = paramD['numThreads']
    
    neighborTL = createNeighborL(geneNames,geneOrderT,synWSize)
    scoresO.initializeScoreArray('synSc') # array to store final synSc result in
    
    ## Prepare argument list
    
    # make list of groups of arguments to be passed to p.imap. There
    # should be numThreads groups.
    argumentL = [[[],tree,neighborTL,numSynToTake,geneNames] for i in range(numThreads)]

    i=0
    for gn1,gn2 in scoresO.iterateEdgesByEndNodes():
        argumentL[i%numThreads][0].append((gn1,gn2))
        i+=1

    ## prepare raw arrays to share
    sharedScoresO.createArrays(scoresO,paramD)

    rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen = sharedScoresO.returnArrays()
    
    ## Run
    with Pool(processes=numThreads,initializer=synScoreGroupInit,initargs=(rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen)) as p:
        for synScoresL in p.imap_unordered(synScoreGroup, argumentL):
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

def synScoreGroupInit(rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen):
    '''Initializer for each separate process doing the synSc
calculation. Loads the global sharedScoresO object with shared
arrays.'''
    sharedScoresO.insertArrays(rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen)

    
def synScoreGroup(argsT):
    '''Given an argument list, including pairs of genes, calculate synteny
    scores. This function is intended to be called by p.map.
    '''

    edgeL,tree,neighborTL,numSynToTake,geneNames = argsT
    
    outL=[]
    for gn1,gn2 in edgeL:
        outL.append(synScore(sharedScoresO,gn1,gn2,tree,neighborTL,numSynToTake,geneNames))
        
    return outL
        
def synScore(sharedScoresO,gn1,gn2,tree,neighborTL,numSynToTake,geneNames):
    '''Given two genes, calculate a synteny score for them. We are given
    the genes, neighborTL, which contains lists of neighbors for each
    gene. For the two sets of neighbors, we find the numSynToTake top
    pairs, and return the average of their scores. The approach is
    greedy. We find the pair with the best score, add it, then remove
    those genes and iterate.
    '''

    sp1 = geneNames.numToStrainNum(gn1)
    sp2 = geneNames.numToStrainNum(gn2)

    L1 = list(neighborTL[gn1])
    L2 = list(neighborTL[gn2])

    topScL= [0] * numSynToTake # min raw score is 0

    for i in range(numSynToTake):
        ind1,ind2,sc = topScore(L1,L2,sharedScoresO)
        if sc == -float('inf'):
            break
        topScL[i] = sc
        del(L1[ind1]) # remove this index
        del(L2[ind2])
        
    synSc = sum(topScL) / numSynToTake
    
    return gn1, gn2, synSc

def topScore(L1,L2,sharedScoresO):
    '''Find the best norm score between genes in L1 and L2. Return the index of
each and the score.'''
    besti1 = 0
    besti2 = 0
    bestSc = -float('inf')

    for i1,gn1 in enumerate(L1):
        for i2,gn2 in enumerate(L2):
            sc = sharedScoresO.getScoreByEndNodes(gn1,gn2,'rawSc')
            if sc != None and sc > bestSc:
                bestSc = sc
                besti1 = i1
                besti2 = i2
    return besti1,besti2,bestSc


#### Calculating the set of all around best reciprocal hit homologs
#    (used in core synteny scores below)

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
    hits1D = getHits(blastDir+strainName1+'_-VS-_'+strainName2+'.out', evalueThresh)
    # ...and 2 blasted against 1
    hits2D = getHits(blastDir+strainName2+'_-VS-_'+strainName1+'.out', evalueThresh)

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


#### Core synteny scores

def calcCoreSynScores(scoresO,strainNamesL,paramD,geneNames,geneOrderT):
    '''Calculate synteny scores based on core genes given in
aabrhL. Scores are between 0 and 1, giving the percentage of syntenic
genes shared.'''

    blastFilePath = paramD['blastFilePath']
    evalueThresh = paramD['evalueThresh']
    aabrhFN = paramD['aabrhFN']
    coreSynWsize = paramD['coreSynWsize']
    
    aabrhL = createAabrhL(blastFilePath,strainNamesL,evalueThresh,aabrhFN)
    
    geneToAabrhT = createGeneToAabrhT(aabrhL,geneNames)
    coreSyntenyT = createCoreSyntenyT(geneToAabrhT,geneOrderT,coreSynWsize)

    # Not parallelized. Pretty fast already.

    # loop over all edges in scoresO, adding corresponding core syn score
    scoresO.initializeScoreArray('coreSynSc') # create array
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


#### Score I/O

def writeScores(scoresO,geneNames,scoresFN):
    '''Write a scores object to file. If scoresFN has the .bout extension, write
binary version, otherwise write in text output format.'''

    scoreTypeL = ['rawSc','synSc','coreSynSc']
    if scoresFN.split('.')[-1] == 'bout':
        scoresO.writeScoresBinary(scoreTypeL,scoresFN)
    else:
        scoresO.writeScoresText(geneNames,scoreTypeL,scoresFN)


def readScores(scoresFN,geneNames=None):
    '''Read scores from file creating a Score object of scores. If
scoresFN has the .bout extension, read binary, otherwise read text
format.
    '''

    scoreTypeL = ['rawSc','synSc','coreSynSc']
    if scoresFN.split('.')[-1] == 'bout':
        scoresO = Score.Score.readScoresBinary(scoreTypeL,scoresFN)
    else:
        scoresO = Score.Score.readScoresText(geneNames,scoreTypeL,scoresFN)
    return scoresO
