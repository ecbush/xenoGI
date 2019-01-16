import sys,glob,os
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import fasta

## funcs

def loadCdsSeq(codingSeqFnL):
    '''Given a list NCBI coding sequence file names (fna) load the sequences and store in a dictionary keyed by locus tag.
    '''
    seqD={}
    for fn in codingSeqFnL:
        for header,seq in fasta.load(fn):
            locusTag = header.split('locus_tag=')[1].split(']')[0]
            seqD[locusTag]=seq
    return seqD

def backAlign(numStrains,protSeqFN, codingSeqD, codeAlignOutFn):
    """prints out blocks of nucleotide alignments based on back alignments from
    the protein alignments."""
    f = open(protSeqFN, 'r')
    S = f.read()
    f.close()
    L = S.split('\n')[:-1] # the way split works, it would produce two
                           # empty strings at end, whereas we're
                           # expecting one below. stip off one newline
    i=0
    outf = open(codeAlignOutFn,'w')
    while i < len(L):
        # do a block
        blockL=[]
        for j in range(numStrains):
            hd = L[i]
            firstPart,locusTag = hd.split(" ")
            xgiGeneName=firstPart[1:] # get rid of >
            i+=1
            protSeq=L[i]
            blockL.append((xgiGeneName,locusTag,protSeq))
            i+=1
        i+=1 # advance past blank line at end of block
        printNucAlign(blockL, codingSeqD,outf)
    outf.close()

def printNucAlign(blockL, codingSeqD,outf):
    '''Prints the nucleotide alignments given by one block of protein
    alignments.'''
    printBlock=''
    for xgiGeneName,locusTag,protSeq in blockL:
        if locusTag not in codingSeqD:
            print(locusTag)
            return
        codingSeq = codingSeqD[locusTag]
        lenProtein = protLength(protSeq)
        # if the protein and sequence don't match in length, throw it out 
        # we have to check two cases because some proteins include the stop
        # codon and some don't
        if lenProtein * 3 != len(codingSeq) and (lenProtein + 1) * 3 != len(codingSeq):
            return
        # printing the header
        printBlock += ">" + xgiGeneName +" " + locusTag + "\n"

        # adding gaps to the nucleotide sequence corresponding to gaps in 
        # the protein sequence and printing the nucleotide sequence
        printBlock += fixSeq(codingSeq, protSeq) + "\n"
    print(printBlock,file=outf)
        
def protLength(prot):
    """returns the length of an aligned protein minus the number of gaps."""
    result =  0
    for x in prot:
        if x != '-':
            result += 1
    return result

def fixSeq(sequence, protAlignment):
    """takes a nucleotide sequence and adds in gaps based on the corresponding
    protein alignment."""
    result = ""
    x = 0
    y = 0

    # while there are still nucleotides in the sequence
    while y != len(protAlignment):

        # if there is a gap in the protein alignment, add gaps to the nucleotide
        # sequence
        if protAlignment[y] == '-':
            result += "---"
            y += 1

        # if there is an amino acid, use the next codon
        else:
            result += sequence[3*x:3*(x+1)].upper()
            x += 1
            y += 1
    return result


if __name__ == "__main__":

    numStrains = int(sys.argv[1]) # will help avoid errors
    protAlignFN = sys.argv[2]
    codeAlignOutFn = sys.argv[3]
    codingSeqFnL= sys.argv[4:]
    codingSeqD = loadCdsSeq(codingSeqFnL)

    backAlign(numStrains,protAlignFN, codingSeqD, codeAlignOutFn)
