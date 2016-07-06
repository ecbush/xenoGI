import sys, fasta

## Funcs

def readGenomeSummary(fn):
    '''Read in a broad genome summary file, return a dictionary keyed by
gene name with start and stop coords.'''
    f = open(fn,"r")
    f.readline()
    D={}
    while True:
        s=f.readline()
        if s=="":
            break
        L=s.split("\t")
        geneName=L[0]
        st=L[4]
        end=L[5]
        D[geneName]=st+"-"+end
    return D

def getGeneName(aaHd):
    '''Takes a header string from aa fasta file, extracts the gene
identifier.'''
    return aaHd.split("|")[1].split()[0]

def writeSimplifiedAA(aaL,genSumD,outFN):
    '''Write out a multi fasta file with the protein name and start and
stop in headers.'''
    f=open(outFN,"w")
    for i in xrange(len(aaL)):
        geneName=getGeneName(aaL[i][0])
        aasq=aaL[i][1]
        stEnd=genSumD[geneName]
        f.write(">"+geneName+" "+stEnd+"\n")
        f.write(aasq+"\n")
    f.close()



if __name__ == "__main__":

    rawProtFileL=[x.rstrip() for x in open("seq/broad/rawProtFiles.txt",'r')]

    #protFN='seq/broad/escherichia_sp._b646_1_proteins.fasta'

    for protFN in rawProtFileL:
        #print protFN
        stem=protFN.split("proteins.fasta")[0]
        genSumFN=stem+"genome_summary_per_gene.txt"

        aaL = fasta.load(protFN)

        genSumD=readGenomeSummary(genSumFN)

        writeSimplifiedAA(aaL,genSumD,stem+".simp.faa")
