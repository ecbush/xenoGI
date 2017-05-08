import sys
import parameters,genbank,blast,trees,genomes,Score,scores


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

    # object for storing scores
    scoresO=Score.Score()
    scoresO.initializeDataAttributes(paramD['blastFilePath'],geneNames)
    
    ## similarity scores
    scoresO = scores.calcRawScores(paramD['fastaFilePath'],paramD['numThreads'],geneNames,paramD['gapOpen'],paramD['gapExtend'],paramD['matrix'],scoresO)

    """
    ## normalized scores
    scoresO,aabrhL,aabrhRawScoreSummmaryD=scores.calcNormScores(tree,strainNum2StrD,paramD['blastFilePath'],paramD['evalueThresh'],scoresO,geneNames,paramD['aabrhFN'])


    ## synteny scores
    scoresO = scores.calcSynScores(scoresO,aabrhRawScoreSummmaryD,geneNames,geneOrderT,paramD['synWSize'],paramD['numSynToTake'],paramD['numThreads'])

    ## core synteny scores
    scoresO = scores.calcCoreSynScores(scoresO,aabrhL,geneNames,geneOrderT,paramD['coreSynWsize'])
    """
    
    # write scores to file
    scores.writeScores(scoresO,geneNames,paramD['scoresFN'])

