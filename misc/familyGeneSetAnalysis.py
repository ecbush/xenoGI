import sys,os
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,genomes,trees,families,scores,islands
from analysis import *

def readGeneLists(fileName):
    '''Read gene lists from fileName. In that file, each set of genes should be on one line, with genes separated by white space. The first element of a line is the name for the gene set, and the rest are the genes.'''
    f=open(fileName,'r')
    geneTL=[]
    while True:
        s=f.readline()
        if s=='':
            break
        geneTL.append(tuple(s.rstrip().split()))
    f.close()
    return geneTL


if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    geneListFN = sys.argv[2]
    
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    
    # get familyL etc.
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    familyL = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    gene2FamD=createGene2FamD(familyL)

    geneTL = readGeneLists(geneListFN)
    # for now just assuming one line, with name first like for syntenyGeneSetAnalysis.py
    
    for gene in geneTL[0][1:]:
        print(str(gene)+"\t"+str(gene2FamD[geneNames.nameToNum(gene)]))
