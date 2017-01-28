import sys
import genbank,blast,trees,genomes,scores,families,groups


if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    ## load data structures we'll use below
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)

    # an object for gene name conversions
    geneNames = genomes.geneNames(params.geneOrderFN,strainStr2NumD,strainNum2StrD)

    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    geneOrderT=genomes.createGeneOrderTs(params.geneOrderFN,geneNames,subtreeL,strainStr2NumD)

     ## read scores
    rawScoresG = scores.readGraph(params.rawScoresFN,geneNames)
    normScoresG = scores.readGraph(params.normScoresFN,geneNames)
    synScoresG = scores.readGraph(params.synScoresFN,geneNames)

    ## make gene families
    familyT = families.families(tree,subtreeL,geneNames,rawScoresG,normScoresG,synScoresG,params.minNormThresh,params.minSynThresh,params.synAdjustThresh,params.synAdjustExtent,params.familyFN,strainNum2StrD)
    
    ## group gene families
    groups.makeGroups(geneOrderT,geneNames,subtreeL,tree,params.groupScoreThreshold,familyT,params.numThreads,strainNum2StrD,params.groupOutFN)
