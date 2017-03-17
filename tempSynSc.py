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
    rawScoresG = scores.readGraph(paramD['rawScoresFN'],geneNames)

    ## normalized scores
    normScoresG = scores.readGraph(paramD['normScoresFN'],geneNames)
    strainNamesL=sorted([strainNum2StrD[leaf] for leaf in trees.leafList(tree)])
    aabrhL = scores.loadOrthos(paramD['aabrhFN'])
    aabrhRawScoreSummmaryD=scores.getAabrhRawScoreSummmaryD(strainNamesL,aabrhL,rawScoresG,geneNames)

    del rawScoresG # don't need anymore, save RAM
    
    ## synteny scores
    synScoresG = scores.createSynScoresGraph(normScoresG,aabrhRawScoreSummmaryD,geneNames,geneOrderT,paramD['synWSize'],paramD['numSynToTake'],paramD['numThreads'],paramD['synScoresFN'])
