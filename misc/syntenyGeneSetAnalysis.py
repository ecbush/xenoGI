import sys,os,statistics,networkx,itertools
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,trees,genomes,scores

# This is currently configured to use the .bout (binary) format of
# scores files. It wouldn't be a very good idea to do with the text
# format, because it would be slow.

def readGeneLists(fileName):
    '''Read gene lists from fileName. In that file, each set of genes should be on one line, with genes separated by white space.'''
    f=open(fileName,'r')
    geneTL=[]
    while True:
        s=f.readline()
        if s=='':
            break
        geneTL.append(tuple(s.rstrip().split()))
    f.close()
    return geneTL
        

def syntenyGeneSet(geneT,synScoresG,geneNames):
    '''Given a tuple of genes in geneT, calculate the average score
between them. If there are some genes that lack a score, we give it
the minimum possible.'''

    subGrNums = (geneNames.nameToNum(name) for name in geneT)
    subgr = synScoresG.subgraph(subGrNums)

    scSum = sum((subgr.get_edge_data(gn1,gn2)['score'] for gn1,gn2 in subgr.edges()))

    # if there was an edge between every node, there would be
    # len(geneT) choose 2.
    maxPossibleNumEdges = len(list(itertools.combinations(geneT,2)))
    actualNumEdges = len(subgr.edges())

    return scSum,maxPossibleNumEdges,actualNumEdges

def printSyntenyGeneList(geneTL,synScoresG,geneNames):
    '''Given a list of gene sets, and a synteny scores graph, print the syntency score for each.'''
    for geneT in geneTL:
        # first element in geneT is name, which we print but don't
        # pass to syntenyGeneSet

        sc,maxPossibleNumEdges,actualNumEdges = syntenyGeneSet(geneT[1:],synScoresG,geneNames)
        print(geneT[0],format(sc,'.4f'),maxPossibleNumEdges,actualNumEdges)
    

    
if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    geneListFN = sys.argv[2]
    
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)
    synScoresG = scores.readGraph(paramD['synScoresFN'],geneNames)
    
    geneTL = readGeneLists(geneListFN)
    
    printSyntenyGeneList(geneTL,synScoresG,geneNames)
