import sys,glob
import fasta

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

def blockReader(f):
    """given a block fasta file, reads in the first block and prints a list
    with elements that are 2 element lists.  The first element of a sublist
    is the name of the gene and the second element is the sequence."""
    blockL = []
    Str = f.readline()

    # loop until we hit the end of a block
    while Str != "\n":
        strain=[]
        if Str=="":
            return []

        # grab the header
        firstPart,locusTag = Str.rstrip().split(" ")
        xgiGeneName=firstPart[1:]
        
        Str = f.readline() # get sequence, which always follows header
        protSeq=Str.rstrip()
        blockL.append((xgiGeneName,locusTag,protSeq))
        Str = f.readline()
        
    return blockL

def backAlign(protSeqFN, codingSeqD):
    """prints out blocks of nucleotide alignments based on back alignments from
    the protein alignments."""
    f = open(protSeqFN, 'r')

    # Going through each block of the protein alignment file
    flag = True
    while flag:
        blockL = blockReader(f)
        if len(blockL) == 0:
            flag = False
        printNucAlign(blockL, codingSeqD)
    f.close()

def printNucAlign(blockL, codingSeqD):
    """Prints the nucleotide alignments given by one block of protein
    alignments."""

    noMissing = True
    printBlock = ""
    for xgiGeneName,locusTag,protSeq in blockL:

        codeSeq = codingSeqD[locusTag]
        lenProtein = protLength(protSeq)

        # if the protein and sequence don't match in length, throw it out 
        # we have to check two cases because some proteins include the stop
        # codon and some don't
        if lenProtein * 3 != len(codeSeq) and \
                (lenProtein + 1) * 3 != len(codeSeq):
            #sys.stderr.write("The protein was not equal to 3 times the length")
            #sys.stderr.write(" of the nucleotide sequence for the gene: \n")
            #sys.stderr.write(result[0] + "\n")
            global blockCounter
            strCounter = str(blockCounter)
            sys.stderr.write(strCounter + "\n")
            noMissing = False
            
        else:
            # printing the header
            printBlock += ">" + xgiGeneName +" " + locusTag + "\n"

            # adding gaps to the nucleotide sequence corresponding to gaps in 
            # the protein sequence and printing the nucleotide sequence
            printBlock += fixSeq(codeSeq, protSeq) + "\n"
    if noMissing:
        print(printBlock)

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

    numStrains = int(sys.argv[1])
    protAlignFN = sys.argv[2]
    codingSeqFnL= sys.argv[3:]
    codingSeqD = loadCdsSeq(codingSeqFnL)

    backAlign(protAlignFN, codingSeqD)
