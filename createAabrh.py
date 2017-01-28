import sys
import genomes,scores

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
    params = __import__(paramFN.replace('.py', ''))

    strainNamesL = loadStrainNames(params.geneOrderFN)
    aabrhL = scores.createAabrhL(params.blastFilePath,strainNamesL,params.evalueThresh,params.aabrhFN)
