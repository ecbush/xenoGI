import sys
fileNamesL=[ln.rstrip() for ln in open(sys.argv[1])]


def locNameFromFasta(fileName):
    '''Read multifasta, return list of (location, gene name ) tuples.'''
    f=open(fileName,"r")
    outL=[]
    while True:
        Str=f.readline()
        if Str=="":
            break
        if Str[0]==">":
            # this is new geneName, put together previous geneName,seq, and move one
            locus = int(Str.split('locus: ')[-1])
            geneName=Str[1:].split()[0]
            outL.append((locus,geneName))
    f.close()
    outL.sort()
    return(outL)    


if __name__ == "__main__":

    for fileName in fileNamesL:
        strainName = fileName.split('/')[-1].split('_aa.fa')[0]
        locNameT = locNameFromFasta(fileName)

        print(strainName + "\t" + " ".join((gn for loc,gn in locNameT)))
