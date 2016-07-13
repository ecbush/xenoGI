# A version of the PhiGs algorithm  http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# Eliot Bush 7/2016
import sys
from tree import *
from genomes import *

## Functions

def createSimilarityGraph(distancesFN,geneName2NumD):
    '''Read distances from the distances file and use to create network
with genes and nodes and edges representing global alignment score
between proteins with significant similarity.'''
    simGrL = [[] for x in geneName2NumD]

    f = open(distancesFN,'r')

    while True:
        s = f.readline()
        if s == '':
            break
        g1,g2,alSc,percSim=s.split('\t')
        alSc = int(alSc)
        percSim = float(percSim)
        simGrL[geneName2NumD[g1]].append((geneName2NumD[g2],alSc,percSim))
        simGrL[geneName2NumD[g2]].append((geneName2NumD[g1],alSc,percSim))
        # could avoid storing alSc twice by putting it in a dict

    # tupleize
    for i in range(len(simGrL)):
        simGrL[i] = tuple(simGrL[i])
    return tuple(simGrL)

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

def createSeedL(tree,leftS,rightS):
    '''Create a list which has the closest match for each gene on the
opposite side of the tree. e.g. if a gene is in tree[1] then we're
looking for the gene in tree[2] with the closest match. For each gene
we get this closest match, put in a list, sort, and return.'''
    

def cluster(nodeOrderL,subtreeL):
    '''Main function.'''
    
    for t,node,lnode,rnode in nodeOrderL:
        if lnode != None:
            # not a tip
            subtree=subtreeL[node]
            leftS,rightS,outgroupS = createLRSets(subtree,geneNum2NameD,geneName2StrainNumD)
            
            

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

    
    #simGrT = createSimilarityGraph(params.distancesFN,geneName2NumD)


    leftS,rightS,outgroupS = createLRSets(tree[1],geneNum2NameD,geneName2StrainNumD)

    nodeOrderL=createNodeProcessOrderList(tree)
