import sys
import genbank,trees,genomes,scores,families,groups,parameters


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
    rawScoresG = scores.readGraph(paramD['rawScoresFN'],geneNames)
    normScoresG = scores.readGraph(paramD['normScoresFN'],geneNames)
    synScoresG = scores.readGraph(paramD['synScoresFN'],geneNames)

    ## make gene families
    familyT = families.families(tree,subtreeL,geneNames,rawScoresG,normScoresG,synScoresG,paramD['minNormThresh'],paramD['minSynThresh'],paramD['synAdjustThresh'],paramD['synAdjustExtent'],paramD['familyFN'],strainNum2StrD)
    
    ## group gene families
    groups.makeGroups(geneOrderT,geneNames,subtreeL,tree,paramD['groupScoreThreshold'],familyT,paramD['numThreads'],strainNum2StrD,paramD['groupOutFN'])
