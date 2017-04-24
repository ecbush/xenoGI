import sys
import parameters,genbank,trees,genomes,scores,families,groups


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
    scoresG = scores.readGraph(paramD['scoresFN'],geneNames)

    ## make gene families
    outputSummaryF = open(paramD['outputSummaryFN'],'w')
    familyT = families.families(tree,subtreeL,geneNames,scoresG,paramD['minNormThresh'],paramD['minCoreSynThresh'],paramD['minSynThresh'],paramD['synAdjustThresh'],paramD['synAdjustExtent'],paramD['familyFN'],strainNum2StrD,outputSummaryF)
    
    ## group gene families
    groups.makeGroups(geneOrderT,geneNames,subtreeL,tree,paramD['proxThreshL'],familyT,paramD['numThreads'],strainNum2StrD,paramD['groupOutFN'],outputSummaryF)
    outputSummaryF.close()
