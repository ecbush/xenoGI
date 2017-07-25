import sys
import parameters,genomes,trees,families,scores,islands
from analysis import *

## Wrapper functions
#  these should be here, as they assume a bunch of global variables.

def printFam(familyNum,fileF=sys.stdout):
    '''This is a wrapper to provide an easy way to print relevant info on
a family. For ease of use, we take only one argument, assuming all the
other required stuff is available at the top level. amilyNum is the
numerical identifier of a family.
    '''

    print("Family error score (count of possibly misassigned genes):",familyL[familyNum].possibleErrorCt,file=fileF)
    
    print(file=fileF)
    print("Matrix of raw similarity scores [0,1] between genes in the family",file=fileF)
    printScoreMatrix(familyNum,subtreeL,familyL,geneNames,scoresO,'rawSc',fileF)
    print(file=fileF)
    print(file=fileF)

    print(file=fileF)
    print("Matrix of normalized similarity scores between genes in the family",file=fileF)
    printScoreMatrix(familyNum,subtreeL,familyL,geneNames,scoresO,'normSc',fileF)
    print(file=fileF)
    print(file=fileF)
    
    print("Matrix of core synteny scores [0,1] between genes in the family",file=fileF)
    printScoreMatrix(familyNum,subtreeL,familyL,geneNames,scoresO,'coreSynSc',fileF)
    print(file=fileF)
    print(file=fileF)

    print("Matrix of synteny scores between genes in the family",file=fileF)
    printScoreMatrix(familyNum,subtreeL,familyL,geneNames,scoresO,'synSc',fileF)
    print(file=fileF)
    print(file=fileF)

    printOutsideFamilyScores(familyNum,subtreeL,familyL,geneNames,scoresO,fileF)
    print(file=fileF)
    print(file=fileF)

def findIsland(searchStr,fileF=sys.stdout):
    '''Print the gene, family and island associated with searchStr. This
is a wrapper that assumes various required objects are present at the
top level.'''
    L=matchFamilyIsland(geneInfoD,geneNames,gene2FamD,fam2IslandD,searchStr)
    for geneName,fam,isl in L:
        print("<gene:"+str(geneName)+">","<family:"+str(fam)+">","<island:"+str(isl)+">",file=fileF)
    

def printIsland(islandNum,synWSize,fileF=sys.stdout):
    '''Print the island and its genomic context in each species. We
    include synWSize/2 genes in either direction beyond the island.
    '''
    printIslandNeighb(islandNum,synWSize,subtreeL,islandByNodeL,familyL,geneOrderT,gene2FamD,fam2IslandD,geneInfoD,geneNames,strainNum2StrD,fileF)

    
def printIslandsAtNode(nodeStr,fileF=sys.stdout):
    '''This is a wrapper to provide an easy way to print all the islands
at a particular node in the tree. For ease of use, we take only a node
number as argument, assuming all the other required stuff is available
at the top level.
    '''
    node = strainStr2NumD[nodeStr]
    vPrintIslands(islandByNodeL[node],subtreeL,familyL,strainNum2StrD,geneNames,fileF)


def printCoreNonCoreByNode(fileF=sys.stdout):
    '''For each node in the focal clade, print the number of core and
non-core families. That is, for that node, we get all families present
in descendent species. We then look at their mrcas. Those with an mrca
below the node in question are non-core, others are core.'''

    focTree = trees.subtree(tree,strainStr2NumD[paramD['rootFocalClade']])
    focInodesL=trees.iNodeList(focTree)
    familyByNodeL=createFamilyByNodeL(geneOrderT,gene2FamD)
    
    rowL=[]
    rowL.append(['Node','Core','Non-Core','Total','% Non-Core'])
    rowL.append(['----','----','--------','-----','----------'])
    for node in focInodesL:
        nonCore,core=coreNonCoreCtAtNode(tree,node,familyByNodeL,familyL)
        rowL.append([strainNum2StrD[node],str(core),str(nonCore),str(core+nonCore),str(format(nonCore/(core+nonCore),".3f"))])
    printTable(rowL,fileF=fileF)
    print(file=fileF)
    return
    
    
if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    
    # load islands
    islandByNodeL=islands.readIslands(paramD['islandOutFN'],tree,strainStr2NumD)
    
    # get familyL etc.
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    
    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])

    familyL = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    gene2FamD=createGene2FamD(familyL)
    fam2IslandD=createFam2IslandD(islandByNodeL)

    # subtree list
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    nodesL=trees.nodeList(tree)


    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    # scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)
    scoresO.createNodeConnectL(geneNames) # make nodeConnectL attribute

    # calc family error scores
    families.calcErrorScores(familyL,scoresO,paramD['minNormThresh'],paramD['minCoreSynThresh'],paramD['minSynThresh'],paramD['famErrorScoreIncrementD'])
    

