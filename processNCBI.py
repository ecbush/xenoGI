import sys,fasta

## Assumes we're working with species which have a single chromosome. We'll produce a file with start and end coords, but no chrom name.


## Funcs

def getProtName(aaHd):
    '''Takes a header string from aa fasta file, extracts the NP protein
identifier.'''
    return aaHd.split("|")[3].split(".")[0]

def getStartEnd(nucHd):
    '''Take header string from nuc fasta file, extracts start and end coords. Assumes there is only one chromosome.'''
    return nucHd.split(":")[1].split()[0]
    
def writeSimplifiedAA(aaL,nucL,outFN):
    '''Write out a multi fasta file with the protein name and start and
stop in headers.'''
    f=open(outFN,"w")
    for i in xrange(len(aaL)):
        protNm=getProtName(aaL[i][0])
        aasq=aaL[i][1]
        stEnd=getStartEnd(nucL[i][0])
        f.write(">"+protNm+" "+stEnd+"\n")
        f.write(aasq+"\n")
    f.close()

            
if __name__ == "__main__":

    aaFileL=[x.rstrip() for x in open(sys.argv[1],'r')]
    nucFileL=[x.rstrip() for x in open(sys.argv[2],'r')]

    for i in range(len(aaFileL)):
        
        aaL = fasta.load(aaFileL[i])
        nucL = fasta.load(nucFileL[i])
        stem=aaFileL[i].split(".faa")[0]
        writeSimplifiedAA(aaL,nucL, stem + ".simp.faa" )
