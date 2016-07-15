# A version of the PhiGs algorithm  http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# Eliot Bush 7/2016
import sys, networkx
from tree import *
from genomes import *


## Functions

def createSimilarityGraph(distancesFN,geneName2NumD):
    '''Read distances from the distances file and use to create network
with genes and nodes and edges representing global alignment score
between proteins with significant similarity.'''

    G=networkx.Graph()
    for i in range(len(geneName2NumD)): G.add_node(i)
    
    f = open(distancesFN,'r')

    while True:
        s = f.readline()
        if s == '':
            break
        g1,g2,alSc,percSim=s.split('\t')
        alSc = int(alSc)
        percSim = float(percSim)

        G.add_edge(geneName2NumD[g1],geneName2NumD[g2],alSc=alSc)

    return G

def nodeDetails(tree,timeOnBrLeadingHere):
    '''Given a tree, return list of (time,node #, left node #, right node
#) for each node. timeOnBrLeadingHere gives the time which should be
added to the time values in the tuple we're making.
    '''
    if tree[1] == ():
        return [(tree[3] + timeOnBrLeadingHere , tree[0], None, None)]
    else:
        l = nodeDetails(tree[1], tree[3] + timeOnBrLeadingHere)
        r = nodeDetails(tree[2], tree[3] + timeOnBrLeadingHere)
        return l + r + [(tree[3] + timeOnBrLeadingHere , tree[0], tree[1][0], tree[2][0])]

        
def createNodeProcessOrderList(tree):
    '''Given a tree, output a list specifying the order of nodes we should
process in the PhiGs algorithm. For each node we give a tuple
(time,node #, left node #, right node #)
    '''
    L=nodeDetails(tree,0)
    L.sort()
    return L

def createLRSets(tree,geneNum2NameD,geneName2StrainNumD):
    '''For every gene in our data, put it into one of three sets. Left,
right, or outgroup. Genes in the left set are found in a species on
the left branch of tree.'''

    leftSpeciesS=set(leafList(tree[1]))
    rightSpeciesS=set(leafList(tree[2]))
    
    leftS=set()
    rightS=set()
    outgroupS=set()
    for geneNum in range(len(geneNum2NameD)): # all genes
        strain=geneName2StrainNumD[geneNum2NameD[geneNum]]
        if strain in leftSpeciesS:
            leftS.add(geneNum)
        elif strain in rightSpeciesS:
            rightS.add(geneNum)
        else:
            outgroupS.add(geneNum)

    return(leftS,rightS,outgroupS)

def closestMatch(gene,rightS,simG):
    '''Find the closest match to gene among the genes in the set geneS in the graph simG. '''

    bestGene=None
    bestEdgeScore = -float('inf')
    for edge in simG.edges_iter(gene):
        if edge[1] in rightS:
            if simG.get_edge_data(*edge)['alSc'] > bestEdgeScore:
                bestEdgeScore =  simG.get_edge_data(*edge)['alSc']
                bestGene = edge[1]
    return bestEdgeScore, gene, bestGene
    
def createSeedL(leftS,rightS,simG):
    '''Create a list which has the closest match for each gene on the
opposite side of the tree. e.g. if a gene is in tree[1] then we're
looking for the gene in tree[2] with the closest match. For each gene
we get this closest match, put in a list, sort, and return.'''

    seedL=[]

    for gene in leftS:
        seedL.append(closestMatch(gene,rightS,simG))
    for gene in rightS:
        seedL.append(closestMatch(gene,leftS,simG))

    seedL.sort(reverse=True)
    return seedL

def getCluster(seed,simG):
    '''Depth first search to get a single cluster.'''

    alThresh,g1,g2 = seed
    
    alreadySearchedS = set()
    notYetSearchedS = set([g1,g2])

    while len(notYetSearchedS) > 0:

        newGenesS=set()
        for gene in notYetSearchedS:
            for edge in simG.edges_iter(gene):
                newGene=edge[1]
                if newGene not in alreadySearchedS and newGene not in newGenesS and simG.get_edge_data(*edge)['alSc'] > alThresh:
                    newGenesS.add(newGene)
            alreadySearchedS.add(gene)
        notYetSearchedS = newGenesS
    return alreadySearchedS
                
    
def family(nodeOrderL,subtreeL):
    '''Main function.'''
    
    for t,node,lnode,rnode in nodeOrderL:
        if lnode != None:
            # not a tip
            subtree=subtreeL[node]
            leftS,rightS,outgroupS = createLRSets(subtree,geneNum2NameD,geneName2StrainNumD)
            seedL = createSeedL(leftS,rightS,simG)
    

## Main

if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    # load data
    tree,strainStr2NumD,strainNum2StrD = readTree(params.treeFN)
   
    geneName2NumD,geneNum2NameD,geneName2StrainNumD = createGeneDs(params.geneOrderFN,strainStr2NumD)

    # subtree list
    subtreeL=createSubtreeL(tree)
    subtreeL.sort()

    
    simG = createSimilarityGraph(params.distancesFN,geneName2NumD)


    nodeOrderL=createNodeProcessOrderList(tree)


    clusterL=[]
    geneUsedL = [False for x in geneName2NumD]
    
    t,node,lnode,rnode = nodeOrderL[0]
    subtree=subtreeL[node]
    leftS,rightS,outgroupS = createLRSets(subtree,geneNum2NameD,geneName2StrainNumD)
    seedL = createSeedL(leftS,rightS,simG)
    
    c=getCluster(seedL[0],simG)
