import sys
import trees,genomes,scores


if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    # load data
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)
   
    geneName2NumD,geneNum2NameD,geneName2StrainNumD = genomes.createGeneDs(params.geneOrderFN,strainStr2NumD)

    # subtree list
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()

    simG = scores.createSimilarityGraph(params.scoresFN,geneName2NumD)

    geneOrderT=genomes.createGeneOrderTs(params.geneOrderFN,geneName2NumD,subtreeL,strainStr2NumD)

    neighborTL = scores.createNeighborL(geneNum2NameD,geneOrderT,params.synWSize)
    
    synScoresG = scores.createSynScoresGraph(simG,neighborTL,params.numSynToTake,params.numThreads)

    scores.writeG(synScoresG,geneNum2NameD,params.synScoresFN)
