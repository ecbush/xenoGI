import sys,glob,networkx
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

    # graph for storing scores
    scoresG=networkx.Graph()
    for geneNum in geneNames.nums: scoresG.add_node(geneNum)
    
    ## similarity scores
    scoresG = scores.calcRawScores(paramD['blastFilePath'],paramD['fastaFilePath'],paramD['numThreads'],geneNames,paramD['gapOpen'],paramD['gapExtend'],paramD['matrix'],scoresG)

    ## normalized scores
    scoresG,aabrhL,aabrhRawScoreSummmaryD=scores.calcNormScores(tree,strainNum2StrD,paramD['blastFilePath'],paramD['evalueThresh'],scoresG,geneNames,paramD['aabrhFN'])

    ## synteny scores
    scoresG = scores.calcSynScores(scoresG,aabrhRawScoreSummmaryD,geneNames,geneOrderT,paramD['synWSize'],paramD['numSynToTake'],paramD['numThreads'])

    ## core synteny scores
    scoresG = scores.calcCoreSynScores(scoresG,aabrhL,geneNames,geneOrderT,paramD['coreSynWsize'])

    # write scores to file
    scores.writeGraph(scoresG,geneNames,paramD['scoresFN'])
    

