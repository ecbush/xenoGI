import sys
import trees,parameters


if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    ## load data structures we'll use below
    treeO,strainStr2NumDO,strainNum2StrDO = trees.readTreeOLD(paramD['treeFN'])

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
