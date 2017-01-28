import sys,glob
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
    rawScoresG = scores.createRawScoresGraph(params.blastFilePath,params.fastaFilePath,params.numThreads,params.rawScoresFN,geneNames,params.gapOpen,params.gapExtend,params.matrix)

    ## normalized scores
    normScoresG,aabrhRawScoreSummmaryD=scores.createNormScoreGraph(tree,strainNum2StrD,params.blastFilePath,params.evalueThresh,rawScoresG,geneNames,params.aabrhFN,params.normScoresFN)    
    
    ## synteny scores
    synScoresG = scores.createSynScoresGraph(normScoresG,aabrhRawScoreSummmaryD,geneNames,geneOrderT,params.synWSize,params.numSynToTake,params.numThreads,params.synScoresFN)
