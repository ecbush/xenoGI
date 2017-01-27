import parasail,networkx,glob,statistics
from multiprocessing import Pool
import genomes,trees

## raw similarity scores

def createRawScoresGraph(blastFilePath,fastaFilePath,numThreads,rawScoresFN,geneNames,gapOpen, gapExtend, matrix):
    '''Find gene pairs with significant blast hits, and get global
alignment scores for each using multiple threads.'''

    # get list of blast files
    blastFnL=glob.glob(blastFilePath)
    protFnL=glob.glob(fastaFilePath)

    # load sequences
    seqD=genomes.loadProt(protFnL)
    
    # Run through blast files getting list of all genes to compare. The
    # blast files have redundancies (ie , g1,g2 and also g2,g1). Keep
    # comparisons to do in set, and don't add if we already have it in
    # either order.
    pairsToDoS = set()
    for fn in blastFnL:
        f = open(fn,'r')
        while True:
            s = f.readline()
            if s=='':
                break
            L = s.split('\t')
            if len(L) != 12: # we only want lines with 12 columns
                continue

            g1 = L[0]
            g2 = L[1]

    
            if (g1,g2) not in pairsToDoS and (g2,g1) not in pairsToDoS:
                # we haven't got it yet, add it
                pairsToDoS.add((g1,g2))
        f.close()

                
    # make list of sets of arguments to be passed to p.map. There
    # should be numThreads sets.
    argumentL = [([],seqD,gapOpen, gapExtend, matrix) for i in range(numThreads)]
    for i,pair in enumerate(pairsToDoS):
        argumentL[i%numThreads][0].append(pair)

        
    # run
    p=Pool(numThreads)
    scoresLL = p.map(rawScoreGroup, argumentL)

    # convert to Graph
    G=networkx.Graph()
    for geneNum in geneNames.nums: G.add_node(geneNum)

    for scoresL in scoresLL:
        for g1,g2,sc in scoresL:
            G.add_edge(geneNames.nameToNum(g1),geneNames.nameToNum(g2),score=sc)

    # write scores to file
    writeGraph(G,geneNames,rawScoresFN)
    
    return G

    
def rawScoreGroup(argT):
    '''Given a dictionary of sequences and a list of gene pairs, go
through each pair and get a needleman wunch based score.
    '''
    pairL,seqD,gapOpen, gapExtend, matrix = argT
    matrix = eval(matrix) # turn it from string to parasail matrix object
    scoresL = []
    for g1,g2 in pairL:
        scaled = rawScore(seqD[g1],seqD[g2],gapOpen, gapExtend, matrix)
        scoresL.append((g1,g2,scaled))
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


def createNormScoreGraph(strainNamesL,blastFilePath,evalueThresh,rawScoresG,geneNames,aabrhFN,normScoresFN):
    '''Given directory of blast output and a graph of raw similarity
scores, calculate normalized similarity scores by comparing each score
with the range of scores in in all around best reciprocal hits in that
pair of strains.'''

    aabrhL = createAabrhL(blastFilePath,strainNamesL,evalueThresh,aabrhFN)

    aabrhRawScoreSummmaryD=getAabrhRawScoreSummmaryD(strainNamesL,aabrhL,rawScoresG,geneNames)
   
    # make norm scores graph
    # get same nodes (genes) as sim graph
    normScoresG=networkx.Graph()
    for node in rawScoresG.nodes_iter(): normScoresG.add_node(node)

    # loop over each edge, normalizing score and putting in normScoresG
    for gn1,gn2 in rawScoresG.edges_iter():
            data=rawScoresG.get_edge_data(gn1,gn2)

            # find mean,std from aabrhRawScoreSummmaryD.
            gnName1 = geneNames.numToName(gn1)
            sp1,restOfName1 = gnName1.split('-')
            gnName2 = geneNames.numToName(gn2)
            sp2,restOfName1 = gnName2.split('-')
            mean,std = aabrhRawScoreSummmaryD[(sp1,sp2)]
            sc = normScore(data['score'],mean,std)
            normScoresG.add_edge(gn1,gn2,score=sc)

    writeGraph(normScoresG,geneNames,normScoresFN)
 
        
    return normScoresG,aabrhRawScoreSummmaryD


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

        # we only want the first hit for any query gene. EB CHANGE THIS LATER to check to get best score even if doesn't come first.
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

def getAabrhRawScoreSummmaryD(strainNamesL,aabrhL,rawScoresG,geneNames):
    '''Given raw scores and a directory with blast output, finds the sets of all around best reciprocal hits. Then for each pair of species, calculates the mean and standard deviation of scores and stores in a dictionary.'''

    # now loop through these, sorting scores into a dict keyed by species pair.

    # create dictionary, (representing an upper triangular matrix)
    spScoreD={}
    for i in range(len(strainNamesL)-1):
        strain1 = strainNamesL[i]
        for j in range(i+1,len(strainNamesL)):
            strain2 = strainNamesL[j]
            spScoreD[(strain1,strain2)]=[]

    # loop through aabrhL and populate
    for orthoT in aabrhL:
        spScoreD = addPairwiseScores(spScoreD,orthoT,rawScoresG,geneNames)

    # get mean and standard deviation
    summaryD = {}
    for sp1,sp2 in spScoreD:
        mean = statistics.mean(spScoreD[(sp1,sp2)])
        std = statistics.stdev(spScoreD[(sp1,sp2)])
        summaryD[(sp1,sp2)] = (mean,std)
        summaryD[(sp2,sp1)] = (mean,std)
        
    return summaryD

def addPairwiseScores(spScoreD,orthoT,rawScoresG,geneNames):
    '''Given a dictionary for storing pairwise scores, and ortholog set in
orthoT, and a network of scores, rawScoresG, pull out all species pairs, and
add score for each in appropriate place ins spScoreD.'''

    for i in range(len(orthoT)-1):
        gene1 = orthoT[i]
        sp1,restOfName1=gene1.split('-')
        geneNum1=geneNames.nameToNum(gene1)
        for j in range(i+1,len(orthoT)):
            gene2 = orthoT[j]
            geneNum2=geneNames.nameToNum(gene2)
            sp2,restOfName1=gene2.split('-')
            data=rawScoresG.get_edge_data(geneNum1,geneNum2)
            sc = data['score']
            key = tuple(sorted([sp1,sp2]))
            spScoreD[key].append(sc)
    return spScoreD

def normScore(rawScore,mean,std):
    '''Given a raw score and mean and std of all raw scores, return a
score normalized by std and centered around zero.'''
    norm = ( rawScore - mean ) / std

    # OLD
    # scale norm to [0,1]
    #normMin = -mean/std
    #normMax = mean/std
    #scaledNorm = ( norm - normMin ) / ( normMax - normMin )
    #return scaledNorm

    return norm


## synteny scores

def createSynScoresGraph(G,aabrhRawScoreSummmaryD,geneNames,geneOrderD,synWSize,numSynToTake,numThreads,synScoresFN):
    '''Create a graph with genes as nodes, and edges representing the
synteny score between two genes. We only bother making synteny scores
for those genes that have an edge in G.
    '''

    neighborTL = createNeighborL(geneNames,geneOrderD,synWSize)


    # prepare argument list for map
    argumentL = []
    for gn1,gn2 in G.edges_iter():
        argumentL.append((gn1,gn2,G,neighborTL,numSynToTake,geneNames,aabrhRawScoreSummmaryD))

            
    p=Pool(numThreads) # num threads
    synScoresL = p.map(synScore, argumentL)

    # make synteny graph
    # get same nodes (genes) as sim graph
    synScoresG=networkx.Graph()
    for node in G.nodes_iter(): synScoresG.add_node(node)

    for gn1,gn2,sc in synScoresL:
        synScoresG.add_edge(gn1,gn2,score=sc)

    writeGraph(synScoresG,geneNames,synScoresFN)
        
    return synScoresG

def createNeighborL(geneNames,geneOrderD,synWSize):
    '''Return a list which specifies the neighbors of each gene. Index of
list corresponds to gene number, and the value located at that index
is a tuple of all genes within a synWSize window. e.g. synWSize 5 means we
go 5 genes in either direction.'''

    neighborTL = [None for x in geneNames.nums]

    for contigT in geneOrderD.values():
        for geneNumT in contigT:
            for i in range(len(geneNumT)):
                end = i + synWSize
                st = i-synWSize if i-synWSize>0 else 0 # st can't be less than 0
                L = list(geneNumT[st:end])
                L.remove(geneNumT[i])
                neighborTL[geneNumT[i]] = tuple(L)

    return neighborTL

def synScore(argsT):
    '''Given two genes, calculate a synteny score for them. We are given
    the genes, neighborTL, which contains lists of neighbors for each
    gene. For the two sets of neighbors, we find the numSynToTake top
    pairs, and return the average of their scores. The approach is
    greedy. We find the pair with the best score, add it, then remove
    those genes and iterate.
    '''
    gn1,gn2,G,neighborTL,numSynToTake,geneNames,aabrhRawScoreSummmaryD = argsT

    # get the min possible score for these two species (this is
    # really for the case of using normalized scores, where it
    # varies by species pair.)
    gnName1 = geneNames.numToName(gn1)
    sp1,restOfName1 = gnName1.split('-')
    gnName2 = geneNames.numToName(gn2)
    sp2,restOfName1 = gnName2.split('-')
    mean,std = aabrhRawScoreSummmaryD[(sp1,sp2)]

    minNormScore = normScore(0,mean,std) # min raw score is 0
    
    L1 = list(neighborTL[gn1])
    L2 = list(neighborTL[gn2])


    topScL= [minNormScore] * numSynToTake

    for i in range(numSynToTake):
        ind1,ind2,sc = topScore(L1,L2,G)
        if sc == -float('inf'):
            break
        topScL[i] = sc

    synSc = sum(topScL) / numSynToTake
    
    return gn1, gn2, synSc

def topScore(L1,L2,G):
    '''Find the best score between genes in L1 and L2. Return the index of
each and the score.'''
    besti1 = 0
    besti2 = 0
    bestSc = -float('inf')

    for i1,gn1 in enumerate(L1):
        for i2,gn2 in enumerate(L2):
            data=G.get_edge_data(gn1,gn2)
            if data != None:
                sc = data['score']
                if sc > bestSc:
                    bestSc = sc
                    besti1 = i1
                    besti2 = i2
    return besti1,besti2,bestSc

## Graph I/O

def readGraph(scoresFN,geneNames):
    '''Read scores from the scores file and use to create network
with genes and nodes and edges representing global alignment score
between proteins with significant similarity.'''

    G=networkx.Graph()
    for geneNum in geneNames.nums: G.add_node(geneNum)
    
    f = open(scoresFN,'r')

    while True:
        s = f.readline()
        if s == '':
            break
        g1,g2,sc=s.split('\t')
        sc = float(sc)
        G.add_edge(geneNames.nameToNum(g1),geneNames.nameToNum(g2),score=sc)
    f.close()
    return G

def writeGraph(G,geneNames,fileName):
    '''Given a graph with genes as nodes, write all edges (pairs of genes)
to file in three columns. Gene 1, gene 2 and score.'''

    f=open(fileName,'w')
    
    for gene1Num,gene2Num in G.edges_iter():
        sc = G.get_edge_data(gene1Num,gene2Num)['score']
        print(geneNames.numToName(gene1Num),geneNames.numToName(gene2Num),format(sc,".6f"),sep='\t',file=f)

    f.close()
