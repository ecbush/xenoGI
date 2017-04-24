import parasail,networkx,glob,statistics
from multiprocessing import set_start_method, Pool
import genomes,trees


#### Global variables

scoresG = None # we'll actually put something in it in functions
               # below. The purpose of doing this has to do with
               # efficient memory use in paralellization of synteny
               # scores.

#### Functions


## raw similarity scores

def calcRawScores(blastFilePath,fastaFilePath,numThreads,geneNames,gapOpen, gapExtend, matrix, scoresG):
    '''Find gene pairs with significant blast hits, and get global
alignment scores for each using multiple threads. Store score in scoresG.'''

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
    p.close()
    p.join()

    
    # store in scoresG
    for scoresL in scoresLL:
        for g1,g2,sc in scoresL:
            scoresG.add_edge(geneNames.nameToNum(g1),geneNames.nameToNum(g2),rawSc=sc)

    return scoresG

    
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


def calcNormScores(tree,strainNum2StrD,blastFilePath,evalueThresh,scoresG,geneNames,aabrhFN):
    '''Given directory of blast output and a graph of raw similarity
scores, calculate normalized similarity scores by comparing each score
with the range of scores in in all around best reciprocal hits in that
pair of strains.'''

    strainNamesL=sorted([strainNum2StrD[leaf] for leaf in trees.leafList(tree)])
    aabrhL = createAabrhL(blastFilePath,strainNamesL,evalueThresh,aabrhFN)

    aabrhRawScoreSummmaryD=getAabrhRawScoreSummmaryD(strainNamesL,aabrhL,scoresG,geneNames)
   
    # loop over each edge in scoresG, normalizing score and saving there
    for gn1,gn2 in scoresG.edges_iter():
            rawSc=scoresG.get_edge_data(gn1,gn2)['rawSc']

            # find mean,std from aabrhRawScoreSummmaryD.
            gnName1 = geneNames.numToName(gn1)
            sp1,restOfName1 = gnName1.split('-')
            gnName2 = geneNames.numToName(gn2)
            sp2,restOfName1 = gnName2.split('-')
            mean,std = aabrhRawScoreSummmaryD[(sp1,sp2)]
            normSc = normScore(rawSc,mean,std)
            scoresG.add_edge(gn1,gn2,normSc=normSc)

    return scoresG,aabrhL,aabrhRawScoreSummmaryD


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

def getAabrhRawScoreSummmaryD(strainNamesL,aabrhL,scoresG,geneNames):
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
        spScoreD = addPairwiseScores(spScoreD,orthoT,scoresG,geneNames)

    # get mean and standard deviation
    summaryD = {}
    for sp1,sp2 in spScoreD:
        mean = statistics.mean(spScoreD[(sp1,sp2)])
        std = statistics.stdev(spScoreD[(sp1,sp2)])
        summaryD[(sp1,sp2)] = (mean,std)
        summaryD[(sp2,sp1)] = (mean,std)
        
    return summaryD

def addPairwiseScores(spScoreD,orthoT,scoresG,geneNames):
    '''Given a dictionary for storing pairwise scores, and ortholog set in
orthoT, and a network of scores, scoresG, pull out all species pairs, and
add score for each in appropriate place in spScoreD.'''

    for i in range(len(orthoT)-1):
        gene1 = orthoT[i]
        sp1,restOfName1=gene1.split('-')
        geneNum1=geneNames.nameToNum(gene1)
        for j in range(i+1,len(orthoT)):
            gene2 = orthoT[j]
            geneNum2=geneNames.nameToNum(gene2)
            sp2,restOfName1=gene2.split('-')
            data=scoresG.get_edge_data(geneNum1,geneNum2)
            sc = data['rawSc']
            key = tuple(sorted([sp1,sp2]))
            spScoreD[key].append(sc)
    return spScoreD

def normScore(rawScore,mean,std):
    '''Given a raw score and mean and std of all raw scores, return a
score normalized by std and centered around zero.'''
    norm = ( rawScore - mean ) / std
    return norm

## synteny scores

def calcSynScores(scoresGArg,aabrhRawScoreSummmaryD,geneNames,geneOrderT,synWSize,numSynToTake,numThreads):
    '''Calculate the synteny score between two genes and add to edge
attributes of scoresG. We only bother making synteny scores for those
genes that have an edge in scoresG.
    '''
    # Put the norm scores graph into a global in this namespace so it
    # will be available in other functions to be called below.
    global scoresG
    scoresG = scoresGArg
    
    neighborTL = createNeighborL(geneNames,geneOrderT,synWSize)

    # prepare argument list for map
    argumentL = []
    for gn1,gn2 in scoresG.edges_iter():
        argumentL.append((gn1,gn2,neighborTL,numSynToTake,geneNames,aabrhRawScoreSummmaryD))

    # find size of chunks map should give to individual processes
    chunkSize = int(len(argumentL) / numThreads)
    
    p=Pool(numThreads) # num threads
    synScoresL = p.map(synScore, argumentL, chunksize = chunkSize)
    p.close()
    p.join()
    
    # add to graph
    for gn1,gn2,sc in synScoresL:
        scoresG.add_edge(gn1,gn2,synSc=sc)

    return scoresG

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

def synScore(argsT):
    '''Given two genes, calculate a synteny score for them. We are given
    the genes, neighborTL, which contains lists of neighbors for each
    gene. For the two sets of neighbors, we find the numSynToTake top
    pairs, and return the average of their scores. The approach is
    greedy. We find the pair with the best score, add it, then remove
    those genes and iterate.
    '''

    global scoresG
    gn1,gn2,neighborTL,numSynToTake,geneNames,aabrhRawScoreSummmaryD = argsT

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
        ind1,ind2,sc = topScore(L1,L2,scoresG)
        if sc == -float('inf'):
            break
        topScL[i] = sc
        del(L1[ind1]) # remove this index
        del(L2[ind2])
        
    synSc = sum(topScL) / numSynToTake
    
    return gn1, gn2, synSc

def topScore(L1,L2,G):
    '''Find the best norm score between genes in L1 and L2. Return the index of
each and the score.'''
    besti1 = 0
    besti2 = 0
    bestSc = -float('inf')

    for i1,gn1 in enumerate(L1):
        for i2,gn2 in enumerate(L2):
            data=G.get_edge_data(gn1,gn2)
            if data != None:
                sc = data['normSc']
                if sc > bestSc:
                    bestSc = sc
                    besti1 = i1
                    besti2 = i2
    return besti1,besti2,bestSc


## Core synteny scores

def calcCoreSynScores(scoresG,aabrhL,geneNames,geneOrderT,coreSynWsize):
    '''Calculate synteny scores based on core genes given in
aabrhL. Scores are between 0 and 1, giving the percentage of syntenic
genes shared.'''

    geneToAabrhT = createGeneToAabrhT(aabrhL,geneNames)
    coreSyntenyT = createCoreSyntenyT(geneToAabrhT,geneOrderT,coreSynWsize)

    # Not parallelized. Seems unecessary to.

    # loop over all edges in scoresG, adding corresponding core syn score
    for gn1,gn2 in scoresG.edges_iter():
        coreSynSc=coreSynScore(coreSyntenyT[gn1],coreSyntenyT[gn2],coreSynWsize)
        scoresG.add_edge(gn1,gn2,coreSynSc=coreSynSc)
        
    return scoresG

def createGeneToAabrhT(aabrhL,geneNames):
    '''Create a tuple where the index corresponds to gene number and the value at that location is the number of the aabrh group to which the gene belongs, or None. The aabrh groups are numbered, simply based on the index where they occur in aabrhL.'''

    geneToAabrhL =[None] * len(geneNames.nums)

    for aabrhNum in range(len(aabrhL)):
        for geneName in aabrhL[aabrhNum]:
            geneNum = geneNames.nameToNum(geneName)
            geneToAabrhL[geneNum] = aabrhNum
    return tuple(geneToAabrhL)

def createCoreSyntenyT(geneToAabrhT,geneOrderT,coreSynWsize):
    '''Create and return core synteny tuple. The index of this corresponds to gene number. The value at that index is a tuple of the aabrh numbers for the syntenic genes within coreSynWsize.'''

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
    
## Graph I/O

def writeGraph(G,geneNames,scoresFN):
    '''Write graph G to file. If scoresFN has the .bout extension, write
binary pickle of G, otherwise write in text output format.'''

    if scoresFN.split('.')[-1] == 'bout':
        writeGraphBinary(G,scoresFN)
    else:
        writeGraphText(G,geneNames,scoresFN)

def writeGraphText(G,geneNames,scoresFN):
    '''Given a graph with genes as nodes, write all edges (pairs of genes)
to a tab delimited text file with a header line.
    '''
    # get header line
    for gene1Num,gene2Num in G.edges_iter():
        dataD = G.get_edge_data(gene1Num,gene2Num)
        scoreTypeL=[]
        for key in dataD:
            scoreTypeL.append(key)
        break
    scoreTypeL.sort()
    
    # open file
    f=open(scoresFN,'w')

    # write header
    f.write("\t".join(['gene1','gene2']+scoreTypeL)+'\n')
    
    for gene1Num,gene2Num in G.edges_iter():
        dataD = G.get_edge_data(gene1Num,gene2Num)
        outStrL=[]
        outStrL.append(geneNames.numToName(gene1Num))
        outStrL.append(geneNames.numToName(gene2Num))
        for scoreType in scoreTypeL:
            outStrL.append(format(dataD[scoreType],".6f"))
    
        f.write("\t".join(outStrL)+'\n')

    f.close()

def writeGraphBinary(G,scoresFN):
    '''Given a graph with genes as nodes, pickle and write to scoresFN.'''
    networkx.write_gpickle(G,scoresFN)

def readGraph(scoresFN,geneNames=None):
    '''Read scores from file creating a networkx graph of scores. If
scoresFN has the .bout extension, read binary pickle of graph,
otherwise read text format.'''
    if scoresFN.split('.')[-1] == 'bout':
        G = readGraphBinary(scoresFN)
    else:
        G = readGraphText(scoresFN,geneNames)
    return G
    
def readGraphText(scoresFN,geneNames):
    '''Read scores from a text file of scores and use to create network
with genes and nodes and edges representing global alignment score
between proteins with significant similarity.
    '''

    G=networkx.Graph()
    for geneNum in geneNames.nums: G.add_node(geneNum)
    
    f = open(scoresFN,'r')

    headerL=f.readline().split()
    scoreTypeL=headerL[2:] # first two cols are gene names
    
    while True:
        s = f.readline()
        if s == '':
            break
        lineL=s.split('\t')
        g1 = lineL[0]
        g2 = lineL[1]
        scoreL=lineL[2:]
        for i,sc in enumerate(scoreL):
            kw = {scoreTypeL[i] : float(sc)}
            G.add_edge(geneNames.nameToNum(g1),geneNames.nameToNum(g2),**kw)
    f.close()
    return G

def readGraphBinary(scoresFN):
    '''Read scores from a binary file, containing a pickled networkx graph
of scores.'''
    G = networkx.read_gpickle(scoresFN)
    return G
