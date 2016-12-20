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

    
    ## similarity scores
    simG = scores.createSimScoresGraph(params.blastFilePath,params.fastaFilePath,params.numThreads,params.scoresFN,geneNames,params.gapOpen,params.gapExtend,params.matrix)
    
    ## synteny scores
    synScoresG = scores.createSynScoresGraph(simG,geneNames,geneOrderT,params.synWSize,params.numSynToTake,params.numThreads,params.synScoresFN)

    
    ## make gene families
    familyT = families.families(tree,subtreeL,geneNames,simG,synScoresG,params.minSynThresh,params.synAdjustThresh,params.synAdjustMaxExtent,params.familyFN,strainNum2StrD)

    
    ## group gene families
    groups.makeGroups(geneOrderT,geneNames,subtreeL,tree,params.groupScoreThreshold,familyT,params.numThreads,strainNum2StrD,params.groupOutFN)
