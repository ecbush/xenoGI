import parasail,statistics,sys
from multiprocessing import set_start_method, Pool
from . import genomes,blast,trees,Score

#### Global variables for use with multithreading

sharedScoresO = Score.sharedScore()

#### raw similarity scores

def calcRawScores(paramD,scoresO):
    '''Get a global alignment based raw score for every edge in scoresO.'''

    numProcesses = paramD['numProcesses']
    gapOpen = paramD['gapOpen']
    gapExtend = paramD['gapExtend']
    matrix = paramD['matrix']
    
    # load sequences
    seqD=genomes.loadSeq(paramD,"_prot.fa")

    # make list of sets of arguments to be passed to p.map. There
    # should be numProcesses sets.
    argumentL = [([],seqD,gapOpen, gapExtend, matrix) for i in range(numProcesses)]

    i=0
    for g1,g2 in scoresO.iterateEdgesByEndNodes():
        edgeNum = scoresO.endNodesToEdge(g1,g2)
        edgeT = edgeNum,g1,g2
        argumentL[i%numProcesses][0].append(edgeT)
        i+=1

    # run in multiple processes
    with Pool(processes=numProcesses) as p:
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

def calcSynScores(scoresO,geneOrderD,paramD):
    '''Calculate the synteny score between two genes and add to edge
attributes of scoresO. We only bother making synteny scores for those
genes that have an edge in scoresO.
    '''

    synWSize = paramD['synWSize']
    numSynToTake = paramD['numSynToTake']
    numProcesses = paramD['numProcesses']
    
    neighborTD = createNeighborD(geneOrderD,synWSize)
    scoresO.initializeScoreArray('synSc') # array to store final synSc result in
    
    ## Prepare argument list
    
    # make list of groups of arguments to be passed to p.imap. There
    # should be numProcesses groups.
    argumentL = [[[],neighborTD,numSynToTake] for i in range(numProcesses)]

    i=0
    for gn1,gn2 in scoresO.iterateEdgesByEndNodes():
        argumentL[i%numProcesses][0].append((gn1,gn2))
        i+=1

    ## prepare raw arrays to share
    sharedScoresO.createArrays(scoresO,paramD)

    rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen = sharedScoresO.returnArrays()
    
    ## Run
    with Pool(processes=numProcesses,initializer=synScoreGroupInit,initargs=(rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen)) as p:
        for synScoresL in p.imap_unordered(synScoreGroup, argumentL):
            for gn1,gn2,sc in synScoresL:
                scoresO.addScoreByEndNodes(gn1,gn2,sc,'synSc')

    return scoresO

def createNeighborD(geneOrderD,synWSize):
    '''Return a dict which specifies the neighbors of each gene. Key
corresponds to gene number, and value is a tuple of all genes within a
synWSize window. e.g. synWSize 10 means we go 5 genes in either
direction.
    '''

    lenInEitherDirec = int(synWSize/2)
    
    neighborTD = {}
    
    for contigT in geneOrderD.values():
        for geneNumT in contigT:
            for i in range(len(geneNumT)):
                end = i + lenInEitherDirec
                st = i-lenInEitherDirec if i-lenInEitherDirec>0 else 0 # st can't be less than 0
                L = list(geneNumT[st:end])
                L.remove(geneNumT[i])
                neighborTD[geneNumT[i]] = tuple(L)

    return neighborTD

def synScoreGroupInit(rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen):
    '''Initializer for each separate process doing the synSc
calculation. Loads the global sharedScoresO object with shared
arrays.'''
    sharedScoresO.insertArrays(rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen)

    
def synScoreGroup(argsT):
    '''Given an argument list, including pairs of genes, calculate synteny
    scores. This function is intended to be called by p.map.
    '''

    edgeL,neighborTD,numSynToTake = argsT
    
    outL=[]
    for gn1,gn2 in edgeL:
        outL.append(synScore(sharedScoresO,gn1,gn2,neighborTD,numSynToTake))
        
    return outL
        
def synScore(sharedScoresO,gn1,gn2,neighborTD,numSynToTake):
    '''Given two genes, calculate a synteny score for them. We are given
    the genes, neighborTD, which contains lists of neighbors for each
    gene. For the two sets of neighbors, we find the numSynToTake top
    pairs, and return the average of their scores. The approach is
    greedy. We find the pair with the best score, add it, then remove
    those genes and iterate.
    '''

    L1 = list(neighborTD[gn1])
    L2 = list(neighborTD[gn2])

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

def createAabrhL(strainNamesL,blastFileJoinStr,blastDir,blastExt,evalueThresh,alignCoverThresh,percIdentThresh,aabrhFN):
    '''Get the sets of all around best reciprocal hits.'''

    rHitsL=getAllReciprocalHits(strainNamesL,blastFileJoinStr,blastDir,blastExt,evalueThresh,alignCoverThresh,percIdentThresh)
    
    # get sets of genes in first species that have a reciprocal best
    # hit in each other species. this is to save time in next step.
    possibleGenesD = getPossibleGenesD(rHitsL)

    # then, out of the possible genes, get those that are all
    # pairwise reciprocal bets hits
    aabrhHardCoreL = getAllAroundBRHL(possibleGenesD, rHitsL)

    # write it to file
    f=open(aabrhFN,'w')
    for orthoT in aabrhHardCoreL:
        f.write("\t".join(map(str,orthoT))+"\n")
    f.close()
        
    return aabrhHardCoreL


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
        orthoL.append(tuple(map(int,L)))
            
    f.close()
    return orthoL

def getAllReciprocalHits(strainNamesL,blastFileJoinStr,blastDir,blastExt,evalueThresh,alignCoverThresh,percIdentThresh):
    '''return an upper-diagonal (N-1)xN matrix where each entry
    [i][j] (j > i) contains a dictionary of best reciprocal hits
    between species i and species j, as indexed in list strainNamesL;
    key is always gene in species i, and value gene in species j'''

    blastD,strainPairL = blast.createBlastD(strainNamesL,blastFileJoinStr,blastDir,blastExt,evalueThresh,alignCoverThresh,percIdentThresh)
    
    rHitsL = []
    
    for i in range(len(strainNamesL)): 
        rHitsL.append([])
        # fill in the rest of the row with the dictionary of best
        # reciprocal hits between species i and j (keyed by species i)
        for j in range(len(strainNamesL)):
            if j > i:
                rHitsL[i].append(getReciprocalHits(blastD[(strainNamesL[i],strainNamesL[j])],blastD[(strainNamesL[j],strainNamesL[i])]))
            else:
                rHitsL[i].append(None)
    return rHitsL

def getReciprocalHits(blastL1,blastL2):
    '''Given matching blast output (in opposite directions) return a
dictionary of reciprocal best hits.'''

    bestHits1D = blast.getBestHitsDictionary(blastL1)
    bestHits2D = blast.getBestHitsDictionary(blastL2)
    
    # then store only the reciprocal best hits
    recipHitsD = {}
    for query in bestHits1D:
        hit = bestHits1D[query][0]
        if hit in bestHits2D:
            nextHit = bestHits2D[hit][0]
            if nextHit == query:
                recipHitsD[query] = hit

    return recipHitsD

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
    
    aabrhHardCoreL = []

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
            aabrhHardCoreL.append((gene,)+hitsT)
    return aabrhHardCoreL


#### Core synteny scores

def calcCoreSynScores(scoresO,strainNamesL,paramD,geneOrderD):
    '''Calculate synteny scores based on core genes given in
aabrhHardCoreL. Scores are between 0 and 1, giving the percentage of syntenic
genes shared.'''

    blastFileJoinStr = paramD['blastFileJoinStr']
    blastFilePath = paramD['blastFilePath']
    blastDir,blastExt = paramD['blastFilePath'].split("*")
    evalueThresh = paramD['evalueThresh']
    alignCoverThresh = paramD['alignCoverThresh']
    percIdentThresh = paramD['percIdentThresh']
    aabrhFN = paramD['aabrhFN']
    coreSynWsize = paramD['coreSynWsize']
    
    aabrhHardCoreL = createAabrhL(strainNamesL,blastFileJoinStr,blastDir,blastExt,evalueThresh,alignCoverThresh,percIdentThresh,aabrhFN)

    geneToAabrhD = createGeneToAabrhD(aabrhHardCoreL)
    coreSyntenyD = createCoreSyntenyD(geneToAabrhD,geneOrderD,coreSynWsize)

    # Not parallelized. Pretty fast already.

    # loop over all edges in scoresO, adding corresponding core syn score
    scoresO.initializeScoreArray('coreSynSc') # create array
    for gn1,gn2 in scoresO.iterateEdgesByEndNodes():
        coreSynSc=coreSynScore(coreSyntenyD[gn1],coreSyntenyD[gn2])
        scoresO.addScoreByEndNodes(gn1,gn2,coreSynSc,'coreSynSc')
        
    return scoresO

def createGeneToAabrhD(aabrhHardCoreL):
    '''Create a dict where key corresponds to gene number and the
value at that location is the number of the aabrh group to which the
gene belongs. The aabrh groups are numbered, simply based on
the index where they occur in aabrhHardCoreL.'''
    geneToAabrhD = {}
    for aabrhNum in range(len(aabrhHardCoreL)):
        for geneNum in aabrhHardCoreL[aabrhNum]:
            geneToAabrhD[geneNum] = aabrhNum
    return geneToAabrhD

def createCoreSyntenyD(geneToAabrhD,geneOrderD,coreSynWsize):
    '''Create and return core synteny tuple. The index of this corresponds
to gene number. The value at that index is a tuple of the aabrh
numbers for the syntenic genes within coreSynWsize.'''

    coreSyntenyD = {}
    for geneNumT,aabrhNumT in iterateGeneAabrhPairsOnContig(geneToAabrhD,geneOrderD):
        for pos in range(len(geneNumT)):
            geneNum = geneNumT[pos]
            coreSyntenyD[geneNum] = getAabrhContext(pos,aabrhNumT,coreSynWsize)
    return coreSyntenyD
                    
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

def iterateGeneAabrhPairsOnContig(geneToAabrhD,geneOrderD):
    '''Iterate over contigs in geneOrderD, yielding matching tuples, one
    with gene numbers, the other with aabrh numbers.
    '''
    for strain,contigT in geneOrderD.items():
        contigL=[]
        for geneNumT in contigT:
            aabrhNumT = convertGeneToAabrhTuple(geneNumT,geneToAabrhD)
            yield geneNumT,aabrhNumT

def convertGeneToAabrhTuple(geneNumT,geneToAabrhD):
    '''Take a tuple of gene numbers and make a corresponding tuple where
every gene number is replaced by its corresponding aabrh number, or
None if it has none.'''
    aabrhNumL=[]
    for geneNum in geneNumT:
        if geneNum in geneToAabrhD:
            aabrhNumL.append(geneToAabrhD[geneNum])
        else:
            # this gene is not part of an aabrh set. put in None
            aabrhNumL.append(None)
    return tuple(aabrhNumL)

def coreSynScore(synT1,synT2):
    '''Given two tuples representing the core genes surrounding each of
two genes, determine the number of shared core genes in each. Return
this, divided by the minimum of the two tuple lengths. Values will be
between 0 and 1.

    '''
    minNumSurrounding = min(len(synT1),len(synT2))
    if minNumSurrounding == 0:
        # avoid div by 0
        return 0
    else:
        ct=0
        for cgene in synT1:
            if cgene in synT2:
                ct+=1
        return ct/minNumSurrounding

#### Plotting scores

def plotScoreHists(paramD):
    """Wrapper to make pdf of histograms of scores."""

    import matplotlib.pyplot as pyplot
    from matplotlib.backends.backend_pdf import PdfPages
    from . import xenoGI
    
    numBins = 80 # num bins in histograms
    strainNamesT = xenoGI.readStrainInfoFN(paramD['strainInfoFN'])
    genesO = genomes.genes(paramD['geneInfoFN'])
    scoresO = readScores(strainNamesT,paramD['scoresFN'])
    aabrhHardCoreL = loadOrthos(paramD['aabrhFN'])
    
    def scoreHists(outFN,scoresO,numBins,scoreType,genesO,aabrhHardCoreL=None):
        '''Read through a scores file, and separate into all pairwise
comparisons. Then plot hist of each.'''

        # currently, this seems to require a display for interactive
        # plots. would be nice to make it run without that...

        pyplot.ioff() # turn off interactive mode
        with PdfPages(outFN) as pdf:
            for strainPair in scoresO.getStrainPairs():
                fig = pyplot.figure()
                scoresL = getScoresStrainPair(scoresO,strainPair,scoreType,genesO,aabrhHardCoreL)
                pyplot.hist(scoresL,bins=numBins, density = True, range = [0,1])
                pyplot.title(strainPair[0]+'-'+strainPair[1])
                pdf.savefig()
                pyplot.close()

    # plot histograms
    for scoreType,outFN in [('rawSc','rawSc.pdf'),('synSc','synSc.pdf'),('coreSynSc','coreSynSc.pdf'),]:
        scoreHists(outFN,scoresO,numBins,scoreType,genesO)

    for scoreType,outFN in [('rawSc','rawScHardCore.pdf'),('synSc','synScHardCore.pdf'),('coreSynSc','coreSynScHardCore.pdf'),]:
        scoreHists(outFN,scoresO,numBins,scoreType,genesO,aabrhHardCoreL)

def getScoresStrainPair(scoresO,strainPair,scoreType,genesO,aabrhHardCoreL):
    '''Get all scores for strainPair. If aabrhHardCoreL is not None, then
only get for pairs in aabrhHardCoreL.'''

    if aabrhHardCoreL==None:
        return list(scoresO.iterateScoreByStrainPair(strainPair,scoreType))
    else:
        scoreL=[]
        for aabrhT in aabrhHardCoreL:
            geneL = getGenesFromStrainT(aabrhT,strainPair,genesO)
            if strainPair[0] == strainPair[1]:
                # if it's a self self, then we only got one gene
                # back. get score of it vs. self.
                gn1=geneL[0]
                gn2=geneL[0]
            else:
                gn1,gn2 = geneL
            scoreL.append(scoresO.getScoreByEndNodes(gn1,gn2,scoreType))
        return scoreL

def getGenesFromStrainT(geneT,strainT,genesO):
    '''Given a tuple of genes (one from each strain) identify and return the two which are from strain pair.'''
    genesInStrainTL=[]
    for geneNum in geneT:
        if genesO.numToStrainName(geneNum) in strainT:
            genesInStrainTL.append(geneNum)
    return genesInStrainTL



#### Score I/O

def writeScores(scoresO,strainNamesT,scoresFN,genesO=None,geneInfoFN=None):
    '''Write a scores object to file. If scoresFN has the .bout extension, write
binary version, otherwise write in text output format.'''

    scoreTypeL = ['rawSc','synSc','coreSynSc']
    if scoresFN.split('.')[-1] == 'bout':
        scoresO.writeScoresBinary(strainNamesT,scoreTypeL,scoresFN)
    else:
        scoresO.writeScoresText(strainNamesT,scoreTypeL,scoresFN,genesO,geneInfoFN)

def readScores(strainNamesT,scoresFN):
    '''Read scores from file creating a Score object of scores. If
scoresFN has the .bout extension, read binary, otherwise read text
format.
    '''

    scoreTypeL = ['rawSc','synSc','coreSynSc']
    if scoresFN.split('.')[-1] == 'bout':
        scoresO = Score.Score.readScoresBinary(strainNamesT,scoreTypeL,scoresFN)
    else:
        scoresO = Score.Score.readScoresText(strainNamesT,scoreTypeL,scoresFN)
    return scoresO
