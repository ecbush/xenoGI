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

    # create group list
    groupByNodeL=groups.createGroupL(familyStrainT,tree)


    # subtree list
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    

    ## cut off the last one, which will always be core families
    #groupByNodeL=groupByNodeL[:-1]

    print("Number of groups per node before merging: ", ' '.join([str(len(x)) for x in groupByNodeL]))
    
    # create score matrix
    print("Creating score matrix.",file=sys.stderr)
    scoreL=groups.createScoreL(groupByNodeL,adjacencyS,subtreeL,familyStrainT)

    # iteratively merge groups
    print("Beginning merging.",file=sys.stderr)
    
    for i in range(len(groupByNodeL)-1):
        print("  Merging",i,file=sys.stderr)
        groups.mergeGroupsAtNode(groupByNodeL,scoreL,adjacencyS,subtreeL,i,params.groupScoreThreshold,familyStrainT)

    print("  Did not merge core groups (last entries in groupByNodeL).")
    print("Merging complete.",file=sys.stderr)

    print("Number of groups per node after merging: ", ' '.join([str(len(x)) for x in groupByNodeL]))

    # write groups
    groups.writeGroups(groupByNodeL,params.groupOutFN)

    print("Groups written.",file=sys.stderr)
