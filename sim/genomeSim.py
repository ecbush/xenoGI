## Genome evolution simulation
## Eliot Bush
import sys,random,copy,glob,os
sys.path.append(os.path.join(sys.path[0],'../'))
import trees, parameters
import pyvolve
from pyvolve import Model, Genetics
import numpy as np
from scipy import linalg


class Genome:
    # for simulating the evolution of a single chromsome
    
    def __init__(self, geneL, seqD, paramD ):
        '''Takes a list of genes (numerical) and a corresponding dict
        of sequences (keyed by number).'''
        self.geneL=geneL
        self.seqD = seqD
        self.paramD = paramD
        
    def evolve(self,branchLen,geneCounter,aaFreq,model):
        '''Evolve along branch for distance branchLen. Return new genome
object and list of events.'''

        geneL = copy.deepcopy(self.geneL)
        seqD = copy.deepcopy(self.seqD)
        logL=[]

        ## pyvolve stuff
        instantaneous_matrix = model.extract_rate_matrix()
        stepBrLen = self.paramD['step']
        transition_matrix = linalg.expm( np.multiply(instantaneous_matrix, stepBrLen) )
        ##
        
        numSteps = int(round(branchLen/paramD['step']))
        for i in range(numSteps):

            # gene del
            geneL,seqD,L = self.deletion(geneL,seqD)
            logL.extend(L)

            # gene dup
            geneL,seqD,L,geneCounter=self.duplication(geneL,seqD,geneCounter)
            logL.extend(L)

            # chr inversion
            geneL,L=self.inversion(geneL)
            logL.extend(L)
            
            # xtrans
            geneL,seqD,L,geneCounter=self.xtrans(geneL,seqD,geneCounter,aaFreq)
            logL.extend(L)

            # seq change
            geneL,seqD = self.mutSeqs(geneL,seqD,transition_matrix,aaFreq)
            
        return Genome(geneL, seqD, self.paramD),geneCounter,logL

    def deletion(self,geneL,seqD):
        '''Delete genes.'''
        logL=[]
        if random.random() < self.paramD['delProb']:
            delSize = random.choice(range(self.paramD['minDel'],self.paramD['maxDel']+1))
            delStart = random.choice(range(0,len(geneL)-delSize))
            delGenes = geneL[delStart:delStart+delSize]
            geneL = geneL[:delStart] + geneL[delStart+delSize:]
            for gene in delGenes:
                del seqD[gene]
            # genes deleted go in second column
            eventStr = 'del\t'+' '.join(map(str,delGenes))
            logL.append(eventStr)
            
        return geneL,seqD,logL

    def duplication(self,geneL,seqD,geneCounter):
        '''Duplicate genes and insert adjacent.'''
        logL=[]
        if random.random() < self.paramD['dupProb']:
            dupSize = random.choice(range(paramD['minDup'],paramD['maxDup']+1))
            dupStart = random.choice(range(0,len(geneL)-dupSize))

            # maybe smaller than dupSize if on end
            genesToCopy = geneL[dupStart:dupStart+dupSize]
            eventStr='dup\t'+' '.join(map(str,genesToCopy))+'\t'
            
            finalGeneL = geneL[:dupStart+dupSize] # start new list of genes

            # add the new genes
            newGenesL=[]
            for i in range(len(genesToCopy)):
                newGenesL.append(geneCounter)
                seqD[geneCounter] = seqD[genesToCopy[i]] # add seq
                geneCounter += 1

            finalGeneL.extend(newGenesL)
            eventStr += ' '.join(map(str,newGenesL))

            finalGeneL.extend(geneL[dupStart+dupSize:])
            geneL = finalGeneL
            logL.append(eventStr)
        return geneL,seqD,logL,geneCounter

    def inversion(self,geneL):
        '''Do a chromosomal inversion.'''
        logL=[]
        if random.random() < self.paramD['invProb']:
            invSize = random.choice(range(paramD['minInv'],paramD['maxInv']+1))
            invStart = random.choice(range(0,len(geneL)-invSize))
            invGenesL = geneL[invStart:invStart+invSize]
            invGenesL.reverse()
            newGeneL = geneL[:invStart] + invGenesL + geneL[invStart+invSize:]
            geneL = newGeneL
            logL.append('inv\t'+' '.join(map(str,invGenesL))) # in new order
        return geneL,logL
        
    def xtrans(self,geneL,seqD,geneCounter,aaFreq):
        '''Horizontal transfer event of unrelated genes from outside the
        clade.
        '''
        logL=[]
        if random.random() < self.paramD['hgtProb']:
            hgtSize = random.choice(range(self.paramD['minHgt'],self.paramD['maxHgt']+1))
            hgtStart = random.choice(range(0,len(geneL)))

            # generate new gene names and put in seqD
            newGenesL=[]
            for i in range(hgtSize):
                newGenesL.append(geneCounter)
                seqLen = random.choice(range(self.paramD['minSeqLen'],self.paramD['maxSeqLen']))
                seqD[geneCounter] = randomSeq(seqLen,aaFreq)
                geneCounter+=1
            finalGeneL = geneL[:hgtStart]
            finalGeneL.extend(newGenesL)
            finalGeneL.extend(geneL[hgtStart:])
            geneL=finalGeneL
            logL.append('hgt\t'+' '.join(map(str,newGenesL)))
        return geneL,seqD,logL,geneCounter

    def mutSeq(self,seq,transition_matrix,aaFreq):
        '''Mutate single sequence acoring to protein eovlution model in
self.paramD['model']. This stuff based on Spielman's suggestions on how
to use pyvolve machinery.
        '''
        # Evolve along lineage
        newIntSeq = np.zeros(len(seq), dtype = int)
        sample_probabilities = np.random.uniform(0, 1, size = len(seq))
        for i in range(len(seq)):
            source = seq[i]
            newIntSeq[i] = sampler( sample_probabilities[i], transition_matrix[source] )

        # consider possibility of small indel
        if random.random() < paramD['smallIndelProb']:
            size = random.choice(range(self.paramD['minIndel'],self.paramD['maxIndel']+1))
            if random.random() < 0.5:
                # del
                start = random.choice(range(len(seq)-self.paramD['maxIndel']))
                newIntSeq = np.concatenate([newIntSeq[:start],newIntSeq[start+size:]])
            else:
                # insertion
                start = random.choice(range(len(seq)))
                insSeq = randomSeq(size,aaFreq)
                newIntSeq = np.concatenate([newIntSeq[:start], insSeq, newIntSeq[start:]])
                
        return newIntSeq

    def mutSeqs(self,geneL,seqD,transition_matrix,aaFreq):
        '''Subject the sequences in seqD to modification according to protein
eovlution model in self.paramD['model'].
        '''
        for gene in geneL:
            seqD[gene] = self.mutSeq(seqD[gene],transition_matrix,aaFreq)
        return geneL,seqD

    def fastaString(self,speciesString):
        '''Return a string in fasta format.'''
        L=[]
        for gene in self.geneL:
            L.append(">"+speciesString+'-'+str(gene))
            L.append(integers_to_protein(self.seqD[gene]))
        return "\n".join(L)
        
def randomResidue(randomValue,freqL,aaL):
    '''Given a randomValue, aaL and a corresponding array of frequencies
figure out which one randomValue corresponds to. In the initial call,
freqL will sum to 1, and randomValue will be between 0 and 1. Returns
a single aa character.
    '''
    if len(aaL) == 1:
        return aaL[0] # catch all, shouldn't be necessary
    elif randomValue < freqL[0]:
        return aaL[0]
    else:
        return randomResidue(randomValue - freqL[0],freqL[1:],aaL[1:])
    
def randomSeq(seqLen,aaFreq):
    '''Generate a random protein sequence seqLen long.'''
    protL=[]
    for i in range(seqLen):
        randomValue = random.random()
        freqL = copy.deepcopy(aaFreq) # make copies becaue radomResidue
        protL.append(randomResidue(randomValue,freqL,AMINO_ACIDS))
    seq= protein_to_integers("".join(protL))
    return seq
    
def randomSeqs(paramD,aaFreq,geneCounter):
    '''Create geneNumber random seqs. Put in seqL. Create corresponding
gene numbers and put in geneL.'''
    geneL = []
    seqL = []
    for i in range(paramD['geneNumber']):
        geneL.append(geneCounter)
        geneCounter += 1
        seqLen = random.choice(range(paramD['minSeqLen'],paramD['maxSeqLen']))
        seqL.append(randomSeq(seqLen,aaFreq))
    return geneL,seqL,geneCounter
        
def createInitialGenome(paramD,aaFreq,geneCounter):
    '''Create an initial genome, geneNumber long. Numbering starts at
value given by geneCounter. Returns a Genome object.'''
    geneL,seqL,geneCounter = randomSeqs(paramD,aaFreq,geneCounter)
    # make seqD
    seqD={}
    for i in range(len(geneL)):
        gene=geneL[i]
        seq=seqL[i]
        seqD[gene]=seq
    return Genome(geneL,seqD,paramD),geneCounter

def sim(tree,genome,paramD,geneCounter,eventLogD,aaFreq,model,strainNum2StrD):
    '''Recursive function to simulate evolution of genome starting from
the root of tree. Returns a list of Genome objects. eventLogD is a
dict (keyed by branch) where we keep track of events in the
simulation.
    '''
    if tree[1] == ():
        tipGenome,geneCounter,logL=genome.evolve(tree[3],geneCounter,aaFreq,model)
        eventLogD[strainNum2StrD[tree[0]]] = logL
        return [(tree[0],tipGenome)],geneCounter
    else:
        # simulate on branch leading to this node
        newGenome,geneCounter,logL=genome.evolve(tree[3],geneCounter,aaFreq,model)
        eventLogD[strainNum2StrD[tree[0]]] = logL

        # left and right subtrees
        tipGenomeL = []
        leftTreeGenomesL,geneCounter = sim(tree[1],newGenome,paramD,geneCounter,eventLogD,aaFreq,model,strainNum2StrD)
        tipGenomeL.extend(leftTreeGenomesL)

        rightTreeGenomesL,geneCounter = sim(tree[2],newGenome,paramD,geneCounter,eventLogD,aaFreq,model,strainNum2StrD)
        tipGenomeL.extend(rightTreeGenomesL)
        
        return tipGenomeL,geneCounter
        

#####
# Some functions to help me use pyvolve, kindly provided by Stephanie Spielman
def protein_to_integers(seq):
    """
        Convert a string of amino acids into an integer list.
    """
    return [ AMINO_ACIDS.index(x) for x in seq ]


def integers_to_protein(ints):
    """
        Convert a list of integers back into amino acids.
    """
    return "".join( [ AMINO_ACIDS[x] for x in ints ] )
    

def sampler(r, probabilities):
    """ 
        Sample an evolved stated based on an array of transition probabilities.
    """
    i = 0
    sum = probabilities[i]
    while sum < r:
        i += 1
        sum += probabilities[i]
    return i
####


def writeLog(tree,strainNum2StrD,logD,fileName):
    '''Write the contents of logD, indicating which branch things come from. The logD has logs for each branch. The entries withing a branch are in order. However we do need to make sure the branches come in a sensible order (pre-order).'''

    # get list of branches. branch gets name of node it leads to
    branchL = [strainNum2StrD[brNum] for brNum in trees.nodeList(tree)] 

    f = open(fileName,'w')
    for branch in branchL:
        L=logD[branch]
        if L != []:
            # skip branches where nothing happened
            L=[str(branch)+'\t'+x for x in L]
            print("\n".join(L),file=f)
    f.close()

def writeGenome(genome,speciesString,fileName):
    '''Given the Genome object genome, write to the file fileName in fasta
format.
    '''
    f = open(fileName,'w')
    print(genome.fastaString(speciesString),file=f)
    f.close()

def writeTipSeqs(tipGenomesL,strainNum2StrD,fastaOutFileDir):
    '''Write the genomes of the final species to file in fasta format.'''

    for branchNum,genome in tipGenomesL:
        branchName = strainNum2StrD[branchNum]
        writeGenome(genome,branchName,fastaOutFileDir+branchName+'.fa')
    
## Main

if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    
    eventLogD = {}
    geneCounter = 0

    
    # get aa freqs to use as stat freq and for starting seqs
    f = pyvolve.EmpiricalModelFrequencies(paramD['model'])
    aaFreq = f.compute_frequencies()
    AMINO_ACIDS = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
    model = Model(paramD['model'])

    
    initialGenome,geneCounter = createInitialGenome(paramD,aaFreq,geneCounter)
    
    tipGenomesL,geneCounter=sim(tree,initialGenome,paramD,geneCounter,eventLogD,aaFreq,model,strainNum2StrD)


    ## Output
    
    # if directory for fastas doesn't exist yet, make it
    fastaDir = paramD['fastaFilePath'].split('*')[0]
    if glob.glob(fastaDir)==[]:
        os.mkdir(fastaDir)
    

    writeLog(tree,strainNum2StrD,eventLogD,paramD['logFile'])
    writeTipSeqs(tipGenomesL,strainNum2StrD,fastaDir)
