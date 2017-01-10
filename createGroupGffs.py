import sys
import trees, families, groups


if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))


    ## load data structures we'll use below
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)

    # load groups
    groupByNodeL=groups.readGroups(params.groupOutFN,tree,strainStr2NumD)
