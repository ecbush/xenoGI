import sys, fasta

## Funcs

def readGenomeSummary(fn):
    '''Read in a broad genome summary file, return a dictionary keyed by
gene name with start and stop coords. Also return a dict keyed by name
with gene description as value.
    '''
    f = open(fn,"r")
    f.readline()
    posD={}
    descripD={}
    while True:
        s=f.readline()
        if s=="":
            break
        L=s.split("\t")
        geneName=L[0]
        st=L[4]
        end=L[5]
        description = L[7]
        posD[geneName]=st+"-"+end
        descripD[geneName]=description
    return posD,descripD

def getGeneName(aaHd):
    '''Takes a header string from aa fasta file, extracts the gene
identifier.'''
    return aaHd.split("|")[1].split()[0]

def writeSimplifiedAA(aaL,genSumD,outFN):
    '''Write out a multi fasta file with the protein name and start and
stop in headers.'''
    f=open(outFN,"w")
    for i in range(len(aaL)):
        geneName=getGeneName(aaL[i][0])
        aasq=aaL[i][1]
        stEnd=genSumD[geneName]
        f.write(">"+geneName+" "+stEnd+"\n")
        f.write(aasq+"\n")
    f.close()



if __name__ == "__main__":

    rawProtFileL=[x.rstrip() for x in open(sys.argv[1],'r')]

    descripF = open(sys.argv[2],'w') 
    
    for protFN in rawProtFileL:
        stem=protFN.split("proteins.fasta")[0]
        genSumFN=stem+"genome_summary_per_gene.txt"

        aaL = fasta.load(protFN)

        posD,descripD=readGenomeSummary(genSumFN)

        writeSimplifiedAA(aaL,posD,stem+".simp.faa")

        for nm,descrip in descripD.items():
            print(nm,descrip,sep="\t",file=descripF)

    descripF.close()
