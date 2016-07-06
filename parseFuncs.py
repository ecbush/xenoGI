# Alex Eng, Bo Lee (??)
import sys

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

        # grab the header and the sequence, getting rid of the newlines
        strain += [Str[:-1]]
        Str = f.readline()
        sequence = ""
        while Str != '' and Str != "\n" and Str[0] != ">":
            sequence += Str[:-1]
            Str = f.readline()
        strain += [sequence]
        blockL += [strain]
    return blockL

def findSequence(f, name):
    """given a file with genes and sequences and the name of a gene, this
    will go through the file and return the sequence corresponding to the
    gene."""

    result = ""
    line = f.readline()
    # look until we find the header or hit the end of the file
    while line != "":
        if line[:len(name) + 1] == ">" + name:
            line = f.readline()[:-1]
            while line != "" and line[0] != ">":
                result += line
                line = f.readline()[:-1]
            return result
        line = f.readline()

def matrixFileReader(f):
    """given a reciprocal matrix file, reads in the upper triangular matrix and
    produces an upper triangular matrix, represented as a list of lists.  Each 
    column and each row are keyed to a species, and each entry in the matrix 
    is a dictionary holding the best reciprocal hits between the two species."""

    matrix = []
    line = f.readline()
    line = f.readline()
    # loop until we hit the EOF
    while line != "":
        row = []

        # loop until we hit EOF or get to the next row
        while line != "" and line[:-1] != "$$$":
            if line[:-1] == "!!!":
                row.append(None)
                line = f.readline()
            else:
                dic = {}
                line = f.readline()
                
                # loop until we hit EOF or get to the next dictionary in the row
                while line != "" and line[0] != ">" and line[:-1] != "!!!"\
                        and line[:-1] != "$$$":
                    genes = line.split("#%#")
                    gene1 = genes[0]
                    gene2 = genes[1][:-1]
                    dic[gene1] = gene2
                    line = f.readline()
                row.append(dic)
        matrix.append(row)
        if line[:-1] == "$$$":
            line = f.readline()
    return matrix

def posHTransFileReader(f):
    """given a file with blocks representing complete groups for gene orthologs,
    returns a list of lists, where each list has 2 elements.  The first
    element is a list of names of species that are in the complete group.  The
    second element is the list of genes that are in the complete group."""

    genesList = []

    line = f.readline()[:-1]

    # loop until EOF
    while line != "":
        try:
            int(line[0])
        except ValueError:
            genesList.append(line.split('\t'))
        line = f.readline()[:-1]
    return genesList

def treeFileReader(f):
    """given a file with a PHYLIP format tree, returns the entire tree as a
    string."""

    tree = ""
    line = f.readline()
    while line != "":
        tree += line[:-1]
        line = f.readline()
    return tree

def getGeneIndex(gene, filePath):
    """given the gene name to look for and the path to the file with the
    informatin about the gene's location in the genome, returns the start
    and stop indices."""

    start = 0
    stop = 0
    f = open(filePath, 'r')
    line = f.readline()

    # loop until EOF
    while line != "":

        # look for the gene name
        if len(line) >= 10:
            if line[:10] == gene[:10]:
                line = line.split('\t')
                start = int(line[4])
                stop = int(line[5])
                break
        line = f.readline()
    f.close()
    return [start, stop]

def getSpeciesMarker(gene):
    """returns the letter combination at the beginning of a gene name that
    indicates which species it comes from."""

    for x in range(len(gene)):
        try:
            int(gene[x])
            return gene[:x]
        except ValueError:
            pass

def getPercentID(f, gene1, gene2):
    """returns a dictionary containg the genes blasted as keys and the blast
    output as a list as the value."""

    line = f.readline()
    found = False

    # loop until EOF
    while line != "":

        # check the line for the names of the two genes
        line = line.split("\t")
        if line[0] == gene1:
            found = True

            # if they both match, return the percent identity value
            if line[1] == gene2:
                return line[2]
        elif found:
            break
        line = f.readline()

def makePIDDic(f):
    """returns a dictionary keyed by a tuple of two gene names and with values
    of percent identities given by blast output."""

    dic = {}
    line = f.readline()

    # loop until EOF
    while line != "":
        line = line.split("\t")

        # grab the best percent id for each pairwise blast
        try:
            oldValue = dic[(line[0], line[1])]
            if float(line[2]) > oldValue:
                dic[(line[0], line[1])] = line[2]
        except KeyError:
            dic[(line[0], line[1])] = line[2]
        line = f.readline()
    return dic
