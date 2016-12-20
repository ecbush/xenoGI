import sys
import groups,trees,genomes

## Main

if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    # load data
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)

    geneNames = genomes.geneNames(params.geneOrderFN,strainStr2NumD,strainNum2StrD)
    
    familyStrainT = groups.createFamilyStrainT(params.familyFN,tree,geneNames,strainStr2NumD)
    
    adjacencyS = genomes.createAdjacencySet(params.geneOrderFN,geneNames)

    # run
    groups.makeGroups(adjacencyS,tree,params.groupScoreThreshold,familyStrainT,params.numThreads,strainNum2StrD,params.groupOutFN)
