# Bo Lee, Alex Eng (??), EB
import sys, parseFuncs, fasta

#nucFileSuffix = "genes.fasta"

#F = 'allDBs.txt'

protSeqFN = sys.argv[1]
broadNucSeqListFN = sys.argv[2]
ncbiProtSeqListFN = sys.argv[3]
ncbiNucSeqListFN = sys.argv[4]



blockCounter = -1

def backAlign(protSeqFN, geneDic):
    """prints out blocks of nucleotide alignments based on back alignments from
    the protein alignments."""
    f = open(protSeqFN, 'r')

    # Going through each block of the protein alignment file
    flag = True
    while flag:
        block = parseFuncs.blockReader(f)
        global blockCounter
        blockCounter += 1
        if len(block) == 0:
            flag = False
        printNucAlign(block, geneDic)
    f.close()
    
def printNucAlign(block, geneDic):
    """prints the nucleotide alignments given by one block of protein
    alignments."""

    noMissing = True
    printBlock = ""
    for x in block:
        #print x
        name = x[0].split(">")[1]
        sequence = geneDic[name]
        lenProtein = protLength(x[1])

        # if the protein and sequence don't match in length, throw it out 
        # we have to check two cases because some proteins include the stop
        # codon and some don't
        if lenProtein * 3 != len(sequence) and \
                (lenProtein + 1) * 3 != len(sequence):
            #sys.stderr.write("The protein was not equal to 3 times the length")
            #sys.stderr.write(" of the nucleotide sequence for the gene: \n")
            #sys.stderr.write(result[0] + "\n")
            global blockCounter
            strCounter = str(blockCounter)
            sys.stderr.write(strCounter + "\n")
            noMissing = False
            
        else:
            # printing the header
            printBlock += x[0] + "\n"

            # adding gaps to the nucleotide sequence corresponding to gaps in 
            # the protein sequence and printing the nucleotide sequence
            printBlock += fixSeq(sequence, x[1]) + "\n"
    if noMissing:
        print printBlock

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
        

def loadNucSeqD(broadNucSeqListFN,ncbiProtSeqListFN,ncbiNucSeqListFN):
    '''Given file names, load all the nucleotide sequences into a dict,
keyed by gene.'''

    geneD={}
    
    # load broad, where gene names are included in header
    broadFileList = [x.rstrip() for x in open(broadNucSeqListFN,"r")]
    for fn in broadFileList:
        L=fasta.load(fn)
        for hd,sq in L:
            geneNm=hd.split()[0][1:] # lop off >
            geneD[geneNm]=sq

    # load ncbi. here we the gene names are not included in the nuc
    # file. so we load the nuc and the protein files, which have genes
    # in the same order. then we get the gene name from the prot file,
    # and the seq from the nuc file
    ncbiProtSeqFileList = [x.rstrip() for x in open(ncbiProtSeqListFN,"r")]
    ncbiNucSeqFileList = [x.rstrip() for x in open(ncbiNucSeqListFN,"r")]
    for i in range(len(ncbiProtSeqFileList)):
        nucL=fasta.load(ncbiNucSeqFileList[i])
        protL=fasta.load(ncbiProtSeqFileList[i])
        for j in range(len(nucL)):
            sq=nucL[j][1]
            geneNm=protL[j][0].split()[0][1:] # lop off >
            geneD[geneNm]=sq

    return geneD


if __name__ == "__main__":
    geneD = loadNucSeqD(broadNucSeqListFN,ncbiProtSeqListFN,ncbiNucSeqListFN)
    backAlign(protSeqFN, geneD)
