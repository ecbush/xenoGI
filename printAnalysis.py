import sys
import parameters,genomes,trees,families,scores,islands
from analysis import *

    
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
    
    # Print out all islands
    islandsOutF = open(paramD['islandsSummaryFN'],'w')
    vPrintAllIslands(islandByNodeL,tree,paramD['rootFocalClade'],subtreeL,familyL,strainStr2NumD,strainNum2StrD,geneNames,islandsOutF)
    islandsOutF.close()

    # Print species files with all the genes, grouped by contig
    printSpeciesContigs(geneOrderT,paramD['genesFNstem'],paramD['genesFNextension'],geneNames,gene2FamD,fam2IslandD,geneInfoD,familyL,strainNum2StrD)
