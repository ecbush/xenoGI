import sys
fileNamesL=[ln.rstrip() for ln in open(sys.argv[1])]


def locNameFromFasta(fileName):
    '''Read multifasta, return list of gene name strings. The order in the
fasta file is the order on the chromosome.'''
    f=open(fileName,"r")
    outL=[]
    while True:
        Str=f.readline()
        if Str=="":
            break
        if Str[0]==">":
            # this is new geneName
            geneName=Str[1:].rstrip()
            outL.append(geneName)
    f.close()
    return(outL)    


if __name__ == "__main__":

    for fileName in fileNamesL:
        strainName = fileName.split('/')[-1].split('.fa')[0]
        geneNameL = locNameFromFasta(fileName)

        print(strainName + "\t" + " ".join(geneNameL))

