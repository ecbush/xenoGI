import sys,os,statistics,itertools
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,trees,genomes,scores

# This is currently configured to use the .bout (binary) format of
# scores files. It wouldn't be a very good idea to do with the text
# format, because it would be slow.

def readGeneLists(fileName):
    '''Read gene lists from fileName. In that file, each set of genes should be on one line, with genes separated by white space. The first element of a line is the name for the gene set, and the rest are the genes.'''
    f=open(fileName,'r')
    geneTL=[]
    while True:
        s=f.readline()
        if s=='':
            break
        geneTL.append(tuple(s.rstrip().split()))
    f.close()
    return geneTL

def getEdges(geneNumsL,scoresO):
    '''Given a list of gene numbers, get a non-redudant list of the edges
that exist between them.'''
    edgeL = []
    for gn in geneNumsL:
        for edge in scoresO.getConnectionsEdge(gn):
            # remove redudancies
            if edge not in edgeL:
                edgeL.append(edge)
    return edgeL
    
def scoreGeneSet(geneT,scoresO,geneNames,scoreType):
    '''Given a tuple of genes in geneT, calculate the average score
between them. If there are some genes that lack a score, we give it
the minimum possible. The kind of score we use is indicated by
scoreType, which corresponds to the type of score in the scores
graph.
    '''

    geneNumsL = [geneNames.nameToNum(name) for name in geneT]
    edgeL = getEdges(geneNumsL,scoresO)
    
    scSum = sum((scoresO.getScoreByEdge(edge,scoreType) for edge in edgeL))

    # if there was an edge between every node, there would be
    # len(geneT) choose 2.
    maxPossibleNumEdges = len(list(itertools.combinations(geneT,2)))
    actualNumEdges = len(edgeL)

    return scSum,maxPossibleNumEdges,actualNumEdges

def printScoreGeneList(geneTL,scoresO,geneNames,scoreType):
    '''Given a list of gene sets, and a scores graph, print the score for
each.'''
    for geneT in geneTL:
        # first element in geneT is name, which we print but don't
        # pass to scoreGeneSet

        sc,maxPossibleNumEdges,actualNumEdges = scoreGeneSet(geneT[1:],scoresO,geneNames,scoreType)
        print(geneT[0],format(sc,'.4f'),maxPossibleNumEdges,actualNumEdges)
    

    
if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    geneListFN = sys.argv[2]
    
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)
    scoresO.createNodeEdgeL(geneNames)

    geneTL = readGeneLists(geneListFN)
    
    printScoreGeneList(geneTL,scoresO,geneNames,'synSc')
    # Key to scoreType:
    # 'rawSc' - raw similarity score
    # 'normSc' - normalized similarity score
    # 'synSc' - synteny score
    # 'coreSynSc' - core synteny score
