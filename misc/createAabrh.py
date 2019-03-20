import sys, os
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import genomes,scores,parameters

def loadStrainNames(geneOrderFN):
    '''Load gene order file, and keep the first thing on each line, which
is strain name.'''
    strainNamesL=[]
    f=open(geneOrderFN)
    while True:
        s=f.readline()
        if s=='':
            break
        strainNamesL.append(s.split()[0])
    strainNamesL.sort()
    return strainNamesL
    
if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)

    strainNamesL = loadStrainNames(paramD['geneOrderFN'])
    aabrhL = scores.createAabrhL(paramD['blastFilePath'],strainNamesL,paramD['evalueThresh'],paramD['aabrhFN'])
