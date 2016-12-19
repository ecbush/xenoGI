import sys
import families,trees,genomes,scores


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

    synScoresG = scores.createSimilarityGraph(params.synScoresFN,geneNames)
    
    nodeOrderL=families.createNodeProcessOrderList(tree)



    familyL = families.families(nodeOrderL,subtreeL,geneNames,simG,synScoresG,params.minSynThresh,params.synAdjustThresh,params.synAdjustMaxExtent)

    families.printFamilies(familyL,geneNames,strainNum2StrD,params.familyFN)

