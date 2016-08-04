# Functions for loading genes and gene order
import trees,fasta

def loadProt(protFnL):
    '''Given a list of file names of simp.faa type fasta files, load the
sequences and store in a dictionary keyed by protein name.
    '''
    seqD={}
    for fn in protFnL:
        for header,seq in fasta.load(fn):
            gn = header.split()[0][1:]
            seqD[gn]=seq
    return seqD

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

def createGeneDescriptionsD(geneDescriptionsFN):
    '''From given file name, create dictionary with gene names as keys and
descriptions as values.'''
    f = open(geneDescriptionsFN,'r')
    D = {}
    while True:
        s = f.readline()
        if s == '':
            break
        s=s.rstrip()
        L=s.split('\t')
        geneName = L[0]
        if len(L)==2:
            geneDescription = L[1]
        else:
            geneDescription = ''
        D[geneName] = geneDescription
        
    f.close()    
    return D
    
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

def createGeneOrderTs(geneOrderFN,geneName2NumD,subtreeL,strainStr2NumD):
    '''Go though gene order file and get orderings into a set of tuples.'''
    f = open(geneOrderFN,'r')
    geneOrderL=[None for x in range(trees.nodeCount(subtreeL[-1]))] # an index for each node
    while True:
        s = f.readline()
        if s == '':
            break
        s=s.rstrip()
        # note, our gene order format has contigs separated by \t, and
        # genes within them separated by a space character.
        L=s.split('\t')
        strain = L[0]
        contigL=[]
        for contig in L[1:]:
            geneNumT=tuple((geneName2NumD[g] for g in contig.split(' ')))
            contigL.append(geneNumT)
            
        geneOrderL[strainStr2NumD[strain]]=tuple(contigL)
    return tuple(geneOrderL)

