# Functions for loading genes and gene order
import sys,glob
from . import fasta
from . import trees

def loadSeq(paramD,fileEnding,genesS=None):
    '''Given paramD and the type of sequence, load the sequences and store
in a dictionary keyed by gene number. fileEnding is a string, either
"_prot.fa" or "_dna.fa". genesS is a set specifying a subset of genes
we want to keep. If missing, we keep all.
    '''
    fileEndingFound = False
    seqD={}
    for fn in glob.glob(paramD['fastaFilePath']):
        if fileEnding in fn:
            fileEndingFound = True # we've seen at least once
            for header,seq in fasta.load(fn):
                gn = int(header.split("_")[0][1:])
                if genesS == None:
                    seqD[gn]=seq
                elif gn in genesS:
                    seqD[gn]=seq

    if not fileEndingFound:
        # we never saw the file ending, the dict is empty.
        raise OSError("There are no fasta files ending in "+fileEnding)    
    return seqD

class genes:
    def __init__(self, geneInfoFN):
        '''genes object. Keeps track of genes present (organized by strain). Also can optionally load a dictionary which allows us to interconvert from numerical to strain represenatations of a gene.'''

        '''
  For later.
 Contains only a single dictionary to convert
one way, either from geneName to number or visa versa. The method
flipDict inverts keys and values in the dictionary.
        '''

        self.geneRangeByStrainD = {}
        self.numGenes = 0
        self.strainGeneRangeT = ()
        self.geneNumToNameD = None
        self.geneInfoD = None
        
        self.initializeGeneRangeByStrainD(geneInfoFN)
        self.initializeStrainGeneRangeT()

    def initializeGeneRangeByStrainD(self,geneInfoFN):
        '''Given a geneInfo file, fill a dict keyed by strain name with values
giving the range of gene numbers in the strain.'''

        f = open(geneInfoFN,'r')

        # handle first line, which gives first strain
        strainName = f.readline().rstrip()[2:]
        enteringNewSpecies = True
        
        while True:
            s = f.readline()
            if s == '':
                break
            elif s[0] == '#':

                # save previous block
                strainRangeEnd = geneNum+1 # make it the number after the last we saw
                self.geneRangeByStrainD[strainName] = (strainRangeStart,strainRangeEnd)

                # get set for next block
                strainName = s.rstrip()[2:] # strainName for next block
                enteringNewSpecies = True
                
            else:
                geneNum,geneName,commonName,locusTag,proteinId,descrip,chrom,start,end,strand=s.rstrip().split('\t')
                geneNum = int(geneNum)
                if enteringNewSpecies == True:
                    strainRangeStart = geneNum
                enteringNewSpecies = False

        f.close()

        # get last strain
        strainRangeEnd = geneNum+1
        self.geneRangeByStrainD[strainName] = (strainRangeStart,strainRangeEnd)
        self.numGenes = strainRangeEnd

    def initializeStrainGeneRangeT(self):
        '''The numToStrainName function works best with a tuple version of the
strain ranges. This method intializes that.'''
        strainRangeL = [(strain,end) for strain,(start,end) in self.geneRangeByStrainD.items()]
        strainRangeL.sort(key=lambda x: x[1]) # sort by end
        self.strainGeneRangeT = tuple(strainRangeL)

    def initializeGeneNumToNameD(self,geneInfoFN,strainNamesL=None):
        '''Given geneInfo file, fill a gene dict keyed by geneNum with value
geneName. If strainNamesL is provided, only load from those strains
present in it. This is an optional part of this data structure, which
isn't loaded by default.
        '''
        self.geneNumToNameD = self.loadDictFromGeneInfoFile(geneInfoFN,strainNamesL,True)

    def initializeGeneInfoD(self,geneInfoFN,strainNamesL=None):
        '''Read gene info from file, returning a dict keyed by gene name with
    information such as description, start position and so on.
        '''
        self.geneInfoD = self.loadDictFromGeneInfoFile(geneInfoFN,strainNamesL,False)

    def loadDictFromGeneInfoFile(self,geneInfoFN,strainNamesL,onlyNames):
        '''Reads in a dictionary from geneInfo file. Keys are geneNum. If
onlyNames is True, then we make the values geneName, otherwise all
other fields.'''

        D = {}
        f = open(geneInfoFN,'r')
        while True:
            s = f.readline()
            if s == '':
                break
            elif s[0] == '#':
                strainName = s.rstrip()[2:]
                if strainNamesL == None or strainName in strainNamesL:
                    keepGenes = True
                else:
                    keepGenes = False
            else:
                geneNum,geneName,commonName,locusTag,proteinId,descrip,chrom,start,end,strand=s.rstrip().split('\t')
                if keepGenes == True:
                    geneNum = int(geneNum)
                    if onlyNames:
                        D[geneNum] = geneName
                    else:
                        D[geneNum] = geneName,commonName,locusTag,proteinId,descrip,chrom,start,end,strand
                        
        f.close()
        return D

    def numToStrainName(self,geneNum):
        '''Given a gene number, return the strain name corresponding. Uses
binary search to find the where geneNum falls in strainGeneRangeT.'''

        st = 0
        end = len(self.strainGeneRangeT)
        
        while True:
            mid = (st + end) // 2
            if self.strainGeneRangeT[mid][1] <= geneNum:
                st = mid + 1
            else:
                # self.strainGeneRangeT[mid][1] > geneNum, so two
                # possibilities, we're in the range, or we're too high
                if self.strainGeneRangeT[mid-1][1] < geneNum:
                    # we're in the range
                    return self.strainGeneRangeT[mid][0]
                else:
                    end = mid
                    if end<=st:
                        return self.strainGeneRangeT[end][0]

    def numToName(self,geneNumber):
        '''Given gene number, return gene name. Assumes self.geneNumToNameD
has been initialized.'''
        return self.geneNumToNameD[geneNumber]

    def iterGenes(self,strainL=None):
        '''Iterate over genes. If strainL is None (or is not given) iterates
over all genes. If strainL is given, iterates over genes from only the
strains present in strainL.'''

        for strainName,endRange in self.strainGeneRangeT:
            # this will give in order of gene number
            if strainL == None or strainName in strainL:
                for geneNum in self.iterGenesStrain(strainName):
                    yield geneNum

    def iterGenesStrain(self,strainName):
        '''Iterates over all genes in strainName.'''
        rangeStart,rangeEnd = self.geneRangeByStrainD[strainName]
        return range(rangeStart,rangeEnd)

    def numToGeneInfo(self,geneNumber):
        '''Given gene number, return other info about gene. Assumes
self.geneInfoD has been initialized.
        '''
        return self.geneInfoD[geneNumber]
    
    def __len__(self):
        return self.numGenes
    
    def __repr__(self):
        return "<genes object with "+str(len(self))+" genes.>"
        
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

def createGeneProximityD(geneOrderD,geneProximityForGroup):
    '''Go though a gene order tuple pulling out pairs of genes that are
within geneProximityForGroup genes of each other. Store in a dict keyed by
gene pair, with value equal to the distance between the two genes,
measured in number of genes.'''
    geneProximityD = {}
    for contigT in geneOrderD.values():
        for geneNumT in contigT:
            for i in range(len(geneNumT)):
                # slide over genes w/win geneProximityForGroup+1
                # wide.
                getProximityInWindow(geneNumT[i:i+(geneProximityForGroup+1)],geneProximityD)
                    
    return geneProximityD

def createGeneOrderD(geneOrderFN,strainNamesL):
    '''Go though gene order file and get orderings into a set of
tuples. Put in a dict, keyed by strain name. We include only the
strains in strainNamesL. If strainNamesL is None, load all.

    '''
    geneOrderD = {}
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
        if strainNamesL == None or strain in strainNamesL:
            # we only load those strains that are in strainNamesL (or all
            # if strainNamesL is None)
            contigL=[]
            for contig in L[1:]:
                geneNumT=tuple((int(g) for g in contig.split(' ')))
                contigL.append(geneNumT)
            geneOrderD[strain] = tuple(contigL)
    f.close()
    return geneOrderD

