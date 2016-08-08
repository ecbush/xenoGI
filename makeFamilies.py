import sys
import families,trees,genomes,scores


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

    synScoresG = scores.createSimilarityGraph(params.synScoresFN,geneName2NumD)
    
    nodeOrderL=families.createNodeProcessOrderList(tree)

    familyL = families.families(nodeOrderL,subtreeL,geneNum2NameD,geneName2StrainNumD,simG,synScoresG,params.minSynThresh,params.synAdjustThresh,params.synAdjustMaxExtent)

    families.printFamilies(familyL,geneNum2NameD,geneName2StrainNumD,params.familyFN)
