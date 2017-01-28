import sys,glob
import genbank,blast,trees,genomes,scores,families,groups


if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    ## load data structures we'll use below
    geneNames = genomes.geneNames(params.geneOrderFN)
    geneOrderD = genomes.createGeneOrderD(params.geneOrderFN,geneNames)
    strainNamesL = sorted(geneOrderD.keys())

    ## similarity scores
    rawScoresG = scores.createRawScoresGraph(params.blastFilePath,params.fastaFilePath,params.numThreads,params.rawScoresFN,geneNames,params.gapOpen,params.gapExtend,params.matrix)
    
    ## normalized scores
    normScoresG,aabrhRawScoreSummmaryD=scores.createNormScoreGraph(strainNamesL,params.blastFilePath,params.evalueThresh,rawScoresG,geneNames,params.aabrhFN,params.normScoresFN)    
    
    ## synteny scores
    synScoresG = scores.createSynScoresGraph(normScoresG,aabrhRawScoreSummmaryD,geneNames,geneOrderD,params.synWSize,params.numSynToTake,params.numThreads,params.synScoresFN)
