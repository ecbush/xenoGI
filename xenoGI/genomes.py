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
    def __init__(self, geneOrderFN=None):
        '''geneNames object to interconvert between the numerical and string
representation of genes. Contains only a single dictionary to convert
one way, either from geneName to number or visa versa. The method
flipDict inverts keys and values in the dictionary.
        '''

        self.geneD = {}
        self.numGenes = 0
        self.strainGeneRangeT = ()
        self.strainNamesO = None
        
        if geneOrderFN != None:
            self.initializeGeneD(geneOrderFN)

    def initializeGeneD(self,geneOrderFN):
        '''Given a geneOrder file, fill a gene dict keyed by geneName with
value geneNumber.'''

        strainGeneRangeL = []
        num=0
        f = open(geneOrderFN,'r')
        while True:
            s = f.readline()
            if s == '':
                break
            L=s.split()
            strainName = L[0]
            for i in range(1,len(L)): 
                geneName=L[i]
                self.geneD[geneName] = num
                num+=1
            endOfRange = num
            # strain ranges will help us get strain from genes
            strainGeneRangeL.append((strainName,endOfRange))
            
        f.close()
        self.numGenes = num
        self.strainGeneRangeT = tuple(strainGeneRangeL)

    def addStrainNamesO(self,strainNamesO):
        '''Add the strainNamesO object in for convenience. We don't do this when initializing because sometimes we may want a geneNames object before we have a tree or strainNamesO.'''
        self.strainNamesO = strainNamesO
        
    def flipDict(self):
        '''Reverse keys and values in geneD.'''
        newD = {}
        while self.geneD != {}:
            key,value = self.geneD.popitem()
            newD[value] = key
        self.geneD = newD
        
    def nameToNum(self,geneName):
        '''Given gene name, return gene number.'''
        if type(geneName) != str:
            raise ValueError("Input of wrong type. Expected a gene in string form.")
        return self.geneD[geneName]

    def numToName(self,geneNumber):
        '''Given gene number, return gene name.'''
        if type(geneNumber) != int:
            raise ValueError("Input of wrong type. Expected a gene in integer form.")
        return self.geneD[geneNumber]

    def iterGeneNums(self):
        return range(0,self.numGenes)

    def copy(self):
        '''Return a copy of self. We will do a deep copy on geneD.'''

        gnO = geneNames()

        # get deep copy of geneD
        newGeneD={}
        for key,value in self.geneD.items():
            newGeneD[key]=value
        gnO.geneD = newGeneD

        # get rest of attributes
        gnO.numGenes = self.numGenes
        gnO.strainGeneRangeT = self.strainGeneRangeT
        gnO.strainNamesO = self.strainNamesO
        
        return gnO

    def numToStrainName(self,geneNum):
        '''Given a gene number, return the strain name corresponding.'''
        return self.numToStrainNameHelper(geneNum,self.strainGeneRangeT)
        
    def numToStrainNameHelper(self,geneNum,strainGeneRangeT):
        '''Recursive helper function to return the strain name corresponding
to geneNum. Checks if geneNum falls in the first range, if not chops
it and recurses.
        '''
        if strainGeneRangeT == ():
            raise IndexError("geneNum falls outside the range of gene numbers for these strains.")
        else:
            strainName,endOfRange = strainGeneRangeT[0]
            if geneNum < endOfRange:
                # its this one
                return strainName
            else:
                return self.numToStrainNameHelper(geneNum,strainGeneRangeT[1:])

    def numToStrainNum(self,geneNum):
        '''Given a gene number, return the strain number corresponding.'''
        return self.strainNamesO.nameToNum(self.numToStrainName(geneNum))
            
    def nameToStrainName(self,geneName):
        '''Given a gene name, return the strain name corresponding.'''
        return self.numToStrainName(self.nameToNum(geneName))

    def nameToStrainNum(self,geneName):
        '''Given a gene name, return the strain number corresponding.'''
        return self.strainNamesO.nameToNum(self.nameToStrainName(geneName))
    
    def __len__(self):
        return len(self.geneD)
    
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

def createGeneOrderTs(geneOrderFN,geneNamesO,subtreeL,strainNamesO):
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
        if not strain in strainNamesO.strToNumD:
            # we only load those things that are in strainNamesO
            continue
        contigL=[]
        for contig in L[1:]:
            geneNumT=tuple((geneNamesO.nameToNum(g) for g in contig.split(' ')))
            contigL.append(geneNumT)
        geneOrderL[strainNamesO.nameToNum(strain)]=tuple(contigL)
    f.close()
    return tuple(geneOrderL)

