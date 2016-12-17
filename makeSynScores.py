import sys
import trees,genomes,scores


if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    # load data
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)

    geneNames = genomes.geneNames(params.geneOrderFN,strainStr2NumD,strainNum2StrD)
    
    # subtree list
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()

    simG = scores.createSimilarityGraph(params.scoresFN,geneNames)

    geneOrderT=genomes.createGeneOrderTs(params.geneOrderFN,geneNames,subtreeL,strainStr2NumD)

    neighborTL = scores.createNeighborL(geneNames,geneOrderT,params.synWSize)
    
    synScoresG = scores.createSynScoresGraph(simG,neighborTL,params.numSynToTake,params.numThreads)

    scores.writeG(synScoresG,geneNames,params.synScoresFN)
