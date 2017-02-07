import sys,glob
import genbank,blast,trees,genomes,scores,families,groups,parameters


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

    
    ## similarity scores
    rawScoresG = scores.createRawScoresGraph(paramD['blastFilePath'],paramD['fastaFilePath'],paramD['numThreads'],paramD['rawScoresFN'],geneNames,paramD['gapOpen'],paramD['gapExtend'],paramD['matrix'])

    ## normalized scores
    normScoresG,aabrhRawScoreSummmaryD=scores.createNormScoreGraph(tree,strainNum2StrD,paramD['blastFilePath'],paramD['evalueThresh'],rawScoresG,geneNames,paramD['aabrhFN'],paramD['normScoresFN'])    
    
    ## synteny scores
    synScoresG = scores.createSynScoresGraph(normScoresG,aabrhRawScoreSummmaryD,geneNames,geneOrderT,paramD['synWSize'],paramD['numSynToTake'],paramD['numThreads'],paramD['synScoresFN'])
