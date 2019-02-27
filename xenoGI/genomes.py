# Functions for loading genes and gene order
import sys
from . import fasta
from . import trees

def loadProt(protFnL):
    '''Given a list of file names of fasta files with the gene name as
header, load the sequences and store in a dictionary keyed by protein
name.
    '''
    seqD={}
    for fn in protFnL:
        for header,seq in fasta.load(fn):
            gn = header.split()[0][1:]
            seqD[gn]=seq
    return seqD

class geneNames:
    def __init__(self, geneOrderFN,strainNameToNumD,strainNumToNameD):
        '''Create a geneName object with lists of genes and methods to interconvert.'''

        self.strainNameToNumD = strainNameToNumD
        self.strainNumToNameD = strainNumToNameD
        
        # get lists of genes
        names=[]
        num=0
        geneNumToStrainNumD={}
        f = open(geneOrderFN,'r')
        while True:
            s = f.readline()
            if s == '':
                break
            L=s.split()
            strain = L[0]
        
            for i in range(1,len(L)): 
                geneName=L[i]
                names.append(geneName)
                geneNumToStrainNumD[num] = strainNameToNumD[strain]

                num+=1

        f.close()

        namesS=set(names)
        if len(namesS) < len(names):
            raise ValueError("There appear to be redudancies in the gene order file.")

        self.geneNumToStrainNumD = geneNumToStrainNumD
        self.names=tuple(names)
        self.nums=tuple(range(len(names))) # these will be the gene numbers for the genes
        
        # create dictionaries for interconverting
        geneNameToNumD={}
        geneNumToNameD={}
        for i in range(len(self.names)):
            geneName=self.names[i]
            num = self.nums[i]
        
            geneNameToNumD[geneName]=num
            geneNumToNameD[num]=geneName

        self.geneNameToNumD = geneNameToNumD
        self.geneNumToNameD = geneNumToNameD


    def nameToNum(self,geneName):
        return self.geneNameToNumD[geneName]

    def numToName(self,geneNumber):
        return self.geneNumToNameD[geneNumber]

    def numToStrainNum(self,geneNumber):
        return self.geneNumToStrainNumD[geneNumber]

    def nameToStrainNum(self,geneName):
        return self.geneNumToStrainNumD[self.nameToNum(geneName)]

    def numToStrainName(self,geneNumber):
        return self.strainNumToNameD[self.geneNumToStrainNumD[geneNumber]]

    def nameToStrainName(self,geneName):
        return self.strainNumToNameD[self.geneNumToStrainNumD[self.nameToNum(geneName)]]

    def isSameStrain(self,geneNum1,geneNum2):
        '''Tests if two genes, given in numerical form, are in the same
strain. Returns boolean.
        '''
        return self.numToStrainNum(geneNum1) == self.numToStrainNum(geneNum2)
        
    def __repr__(self):
        return "<geneName object with "+str(len(self.names))+" genes.>"
    
def readGeneInfoD(geneInfoFN):
    '''Read gene info from file, returning a dict keyed by gene name with
information such as description, start position and so on.
    '''
    geneInfoD = {}
    f = open(geneInfoFN,'r')
    while True:
        s = f.readline()
        if s == '':
            break
        geneName,commonName,locusTag,descrip,chrom,start,end,strand=s.rstrip().split('\t')
        geneInfoD[geneName]=(commonName,locusTag,descrip,chrom,start,end,strand)
    f.close()
    return geneInfoD

def getProximityInWindow(geneWinT,geneProximityD):
    '''Given a window of genes, calculate the distances between first gene
and the rest (in number of genes) and store in geneProximityD. Updates
geneProximityD in place, returning None.
    '''
    gnA=geneWinT[0]
    for i in range(1,len(geneWinT)):
        gnB=geneWinT[i]

        if gnA<gnB: # always put lower gene number first
            geneProximityD[(gnA,gnB)]=i
        else:
            geneProximityD[(gnB,gnA)]=i

def createGeneProximityD(geneOrderT,geneProximityForGroup):
    '''Go though a gene order tuple pulling out pairs of genes that are
within geneProximityForGroup genes of each other. Store in a dict keyed by
gene pair, with value equal to the distance between the two genes,
measured in number of genes.'''
    geneProximityD = {}
    for contigT in geneOrderT:
        if contigT != None:
            # internal nodes are None, having no genes to be adjacent
            for geneNumT in contigT:
                for i in range(len(geneNumT)):
                    # slide over genes w/win geneProximityForGroup+1
                    # wide.
                    getProximityInWindow(geneNumT[i:i+(geneProximityForGroup+1)],geneProximityD)
                    
    return geneProximityD

def createGeneOrderTs(geneOrderFN,geneNames,subtreeL,strainStr2NumD):
    '''Go though gene order file and get orderings into a set of
tuples. Returns a tuple whose index is strain number, and the value at
that index is a set of tuples representing the contigs.'''
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
            geneNumT=tuple((geneNames.nameToNum(g) for g in contig.split(' ')))
            contigL.append(geneNumT)
            
        geneOrderL[strainStr2NumD[strain]]=tuple(contigL)
    f.close()
    return tuple(geneOrderL)

