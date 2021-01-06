import sys,os
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import parameters,xenoGI


if __name__ == "__main__":

    paramFN = sys.argv[1]
    speciesTreeFN = sys.argv[2] 

    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)
    speciesRtreeO,subtreeD=xenoGI.loadTreeRelatedData(speciesTreeFN)

    strainNamesT = tuple((leaf for leaf in speciesRtreeO.leaves()))
    
    for strainNum in range(len(strainNamesT)):
        print(str(strainNum)+"\t"+strainNamesT[strainNum])
