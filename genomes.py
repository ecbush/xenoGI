# Functions for loading genes and gene order

def createGeneDs(geneOrderFN,strainStr2NumD):
    '''Load the gene order file, and give each gene a unique
number. Returns 3 dictionaries for interconverting between names and
numbers, and for giving strain number given a gene string.
    '''
    num=0
    geneName2NumD={}
    geneNum2NameD={}
    geneName2StrainNumD={}
    f = open(geneOrderFN,'r')
    while True:
        s = f.readline()
        if s == '':
            break
        L=s.split()
        strain = L[0]

        for geneName in L[1:]:
            
            geneName2NumD[geneName]=num
            geneNum2NameD[num]=geneName
            geneName2StrainNumD[geneName]=strainStr2NumD[strain]
            
            num+=1

    f.close()
            
    return geneName2NumD,geneNum2NameD,geneName2StrainNumD

def createAdjacencySet(geneOrderFN,geneName2NumD):
    '''Go though gene order file pulling out pairs of adjacent gene and
putting them in a set.'''
    adjacencyS=set()
    f = open(geneOrderFN,'r')
    while True:
        s = f.readline()
        if s == '':
            break
        s=s.rstrip()
        # note, our gene order format has contigs separated by \t, and
        # genes within them separated by a space character.
        L=s.split('\t')
        strain = L[0]
        for contig in L[1:]:
            geneNameL=contig.split(' ')
            for i in range(len(geneNameL)-1):
                gnA=geneName2NumD[geneNameL[i]]
                gnB=geneName2NumD[geneNameL[i+1]]
                if gnA<gnB: # always put lower gene number first
                    adjacencyS.add((gnA,gnB))
                else:
                    adjacencyS.add((gnB,gnA))
    return adjacencyS

