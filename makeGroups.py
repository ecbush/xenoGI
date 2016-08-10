import sys
import groups,trees,genomes

## Main

if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    # load data
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)

    geneName2NumD,geneNum2NameD,geneName2StrainNumD = genomes.createGeneDs(params.geneOrderFN,strainStr2NumD)

    
    familyStrainT = groups.createFamilyStrainT(params.familyFN,tree,geneName2NumD,geneName2StrainNumD)

    adjacencyS = genomes.createAdjacencySet(params.geneOrderFN,geneName2NumD)


    # run
    groups.makeGroups(adjacencyS,tree,params.groupScoreThreshold,familyStrainT,params.numThreads,params.groupOutFN)
