import sys
import parameters,genbank,trees,genomes,scores,families,islands


if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    ## load data structures we'll use below
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])

    # an object for gene name conversions
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    ## read scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)

    ## make gene families
    outputSummaryF = open(paramD['outputSummaryFN'],'w')
    familyT = families.families(tree,subtreeL,geneNames,scoresO,paramD['minNormThresh'],paramD['minCoreSynThresh'],paramD['minSynThresh'],paramD['synAdjustThresh'],paramD['synAdjustExtent'],paramD['familyFN'],strainNum2StrD,outputSummaryF)
    
    ## group gene families into islands
    islands.makeIslands(geneOrderT,geneNames,subtreeL,tree,paramD['proxThreshL'],familyT,paramD['numThreads'],strainStr2NumD,strainNum2StrD,paramD['rootFocalClade'],paramD['islandOutFN'],outputSummaryF)
    outputSummaryF.close()
