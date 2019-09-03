# Functions for a modified version of the PhiGs algorithm
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# (We've added the use of synteny)
import sys,numpy,random
from scipy.signal import find_peaks
from . import trees,scores
from .Family import *
from .analysis import printTable

#### Main function

def createFamiliesO(tree,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,subtreeD,outputSummaryF):
    '''Given a scoresO object, create gene families using the DTLOR approach.
    '''



## Supporting functions




    
    
## Input/output

def writeFamilies(familiesO,genesO,strainNamesT,paramD):
    '''Write all gene families to fileName, one family per line.'''

    familyFN = paramD['familyFN']
    geneInfoFN = paramD['geneInfoFN']

    # get the num to name dict, only for strains we're looking at.
    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT))
    
    f=open(familyFN,'w')
    for fam in familiesO.iterFamilies():
        f.write(fam.fileStr(genesO)+'\n')
    f.close()


def readFamilies(familyFN,tree,genesO):
    '''Read the family file named familyFN, creating a Families object.
    '''
    familiesO = Families(tree)
    f=open(familyFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.split('\t')
        famNum=int(L[0])
        mrca = L[1]
        if L[2] == "-":
            seedPairL = None
        else:
            seedPairL = [L[2],L[3]]

        lfL = L[4:]

        familiesO.initializeFamily(famNum,mrca,seedPairL)

        for lfStr in lfL:
            lfSplitL = lfStr.rstrip().split(',')
            locusFamNum=int(lfSplitL[0])
            lfMrca = lfSplitL[1]
            geneL=[]
            for geneName in lfSplitL[2:]:
                geneNum = int(geneName.split('_')[0])
                geneL.append(geneNum)
            lfO = LocusFamily(famNum,locusFamNum,lfMrca)
            lfO.addGenes(geneL,genesO)
            familiesO.addLocusFamily(lfO)

    f.close()
    return familiesO
