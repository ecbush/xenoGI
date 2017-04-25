import sys
import genomes,trees,families,scores,groups,parameters
from analysis import *

## Wrapper functions
#  these should be here, as they assume a bunch of global variables.

def printFam(synWSize,family):
    '''This is a wrapper to provide an easy way to print relevant info on
a family. For ease of use, we take only two arguments, assuming all
the other required stuff is available at the top level. Family is the
numerical identifier of a family. synWSize is the size of the genomic
window we will include when printing out genomic context in each
species.
    '''

    print()
    print("Matrix of raw similarity scores [0,1] between genes in the family")
    printScoreMatrix(family,subtreeL,familyT,geneNames,scoresG,'rawSc')
    print()
    print()

    print()
    print("Matrix of normalized similarity scores between genes in the family")
    printScoreMatrix(family,subtreeL,familyT,geneNames,scoresG,'normSc')
    print()
    print()
    
    print("Matrix of core synteny scores between genes in the family")
    printScoreMatrix(family,subtreeL,familyT,geneNames,scoresG,'coreSynSc')
    print()
    print()

    print("Matrix of synteny scores between genes in the family")
    printScoreMatrix(family,subtreeL,familyT,geneNames,scoresG,'synSc')
    print()
    print()

    printOutsideFamilyScores(family,subtreeL,familyT,geneNames,scoresG)
    print()
    print()

    print("Synteny information")
    printFamNeighb(family,synWSize,subtreeL,familyT,geneOrderT,gene2FamD,fam2GroupD,geneInfoD,geneNames,strainNum2StrD)

def printGroupsAtNode(nodeStr):
    '''This is a wrapper to provide an easy way to print all the groups at
a particular node in the tree. For ease of use, we take only a node
number as argument, assuming all the other required stuff is available
at the top level.
    '''
    node = strainStr2NumD[nodeStr]
    vPrintGroups(groupByNodeL[node],subtreeL,familyT,strainNum2StrD,geneNames)

    



if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    
    # load groups
    groupByNodeL=groups.readGroups(paramD['groupOutFN'],tree,strainStr2NumD)
    
    # get familyT etc.
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    
    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])

    familyT = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    
    gene2FamD=createGene2FamD(familyT)
    fam2GroupD=createFam2GroupD(groupByNodeL)

    # subtree list
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()


    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    # scores
    scoresG = scores.readGraph(paramD['scoresFN'],geneNames)
    
