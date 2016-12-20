import parasail,networkx,glob
from multiprocessing import Pool
import genomes

## raw similarity scores

def createSimScoresGraph(blastFilePath,fastaFilePath,numThreads,scoresFN,geneNames,gapOpen, gapExtend, matrix):
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
    scoresLL = p.map(simScoreGroup, argumentL)

    # convert to Graph
    G=networkx.Graph()
    for geneNum in geneNames.nums: G.add_node(geneNum)

    for scoresL in scoresLL:
        for g1,g2,sc in scoresL:
            G.add_edge(geneNames.nameToNum(g1),geneNames.nameToNum(g2),score=sc)

    # write scores to file
    writeGraph(G,geneNames,scoresFN)
    
    return G

    
def simScoreGroup(argT):
    '''Given a dictionary of sequences and a list of gene pairs, go
through each pair and get a needleman wunch based score.
    '''
    pairL,seqD,gapOpen, gapExtend, matrix = argT
    matrix = eval(matrix) # turn it from string to parasail matrix object
    scoresL = []
    for g1,g2 in pairL:
        scaled = simScore(seqD[g1],seqD[g2],gapOpen, gapExtend, matrix)
        scoresL.append((g1,g2,scaled))
    return scoresL
      
   
def simScore(s1,s2,gapOpen, gapExtend, matrix):
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


## synteny scores

def createNeighborL(geneNames,geneOrderT,synWSize):
    '''Return a list which specifies the neighbors of each gene. Index of
list corresponds to gene number, and the value located at that index
is a tuple of all genes within a synWSize window. e.g. synWSize 5 means we
go 5 genes in either direction.'''

    neighborTL = [None for x in geneNames.nums]

    for contigT in geneOrderT:
        if not contigT == None:
            for geneNumT in contigT:
                for i in range(len(geneNumT)):
                    end = i + synWSize
                    st = i-synWSize if i-synWSize>0 else 0 # st can't be less than 0
                    L = list(geneNumT[st:end])
                    L.remove(geneNumT[i])
                    neighborTL[geneNumT[i]] = tuple(L)

    return neighborTL

def pairScore(gn1,gn2,G):
    '''Given a pair of genes, see if there is an edge. If so, return
score, if not, return 0.'''
    data=G.get_edge_data(gn1,gn2)
    if data == None:
        return 0
    else:
        return data['score']


def topScore(L1,L2,simG):
    '''Find the best score between genes in L1 and L2. Return the index of
each and the score.'''
    besti1 = 0
    besti2 = 0
    bestSc = 0
    for i1,gn1 in enumerate(L1):
        for i2,gn2 in enumerate(L2):
            sc = pairScore(gn1,gn2,simG)
            if sc > bestSc:
                bestSc = sc
                besti1 = i1
                besti2 = i2
    return besti1,besti2,bestSc

def synScore(argsT):
    '''Given two genes, calculate a synteny score for them. We are given
    the genes, neighborTL, which contains lists of neighbors for each
    gene. For the two sets of neighbors, we find the numSynToTake top
    pairs, and return their average of their scores. This score is
    between 0 and 1. The approach is greedy. We find the pair with the
    best score, add it, then remove those genes and iterate.
    '''
    gn1,gn2,simG,neighborTL,numSynToTake = argsT
    
    L1 = list(neighborTL[gn1])
    L2 = list(neighborTL[gn2])

    scSum = 0
    counter = 0
    while counter < numSynToTake and len(L1) > 0 and len(L2) > 0:
        ind1,ind2,sc = topScore(L1,L2,simG)
        scSum += sc
        del L1[ind1]
        del L2[ind2]
        counter += 1

    return gn1, gn2, scSum / numSynToTake


def createSynScoresGraph(simG,geneNames,geneOrderT,synWSize,numSynToTake,numThreads,synScoresFN):
    '''Create a graph with genes as nodes, and edges representing the
synteny score between two genes. We only bother making synteny scores
for those genes that have an edge in simG.
    '''

    neighborTL = createNeighborL(geneNames,geneOrderT,synWSize)

    
    p=Pool(numThreads) # num threads
    
    argumentL = []
    # create and edge for every one in simG. Give it weight of sum of
    # scores of adjacent genes
    for gn1,gn2 in simG.edges_iter():
        argumentL.append((gn1,gn2,simG,neighborTL,numSynToTake))

    synScoresL = p.map(synScore, argumentL)

    # make synteny graph
    # get same nodes (genes) as sim graph
    synScoresG=networkx.Graph()
    for node in simG.nodes_iter(): synScoresG.add_node(node)

    for gn1,gn2,sc in synScoresL:
        synScoresG.add_edge(gn1,gn2,score=sc)

    writeGraph(synScoresG,geneNames,synScoresFN)
        
    return synScoresG


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
