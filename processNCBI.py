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

    aaL = fasta.load("seq/ncbi/Escherichia_coli_CFT073.faa")
    nucL = fasta.load("seq/ncbi/Escherichia_coli_CFT073.ffn")
    writeSimplifiedAA(aaL,nucL,"seq/ncbi/Escherichia_coli_CFT073.simp.faa")

    aaL = fasta.load("seq/ncbi/escherichia_coli_O157_H7.faa")
    nucL = fasta.load("seq/ncbi/escherichia_coli_O157_H7.ffn")
    writeSimplifiedAA(aaL,nucL,"seq/ncbi/escherichia_coli_O157_H7.simp.faa")

    aaL = fasta.load("seq/ncbi/Klebsiella_pneumoniae_MGH_78578.faa")
    nucL = fasta.load("seq/ncbi/Klebsiella_pneumoniae_MGH_78578.ffn")
    writeSimplifiedAA(aaL,nucL,"seq/ncbi/Klebsiella_pneumoniae_MGH_78578.simp.faa")

    aaL = fasta.load("seq/ncbi/Salmonella_Typhimurium_14028S.faa")
    nucL = fasta.load("seq/ncbi/Salmonella_Typhimurium_14028S.ffn")
    writeSimplifiedAA(aaL,nucL,"seq/ncbi/Salmonella_Typhimurium_14028S.simp.faa")

    aaL = fasta.load("seq/ncbi/Shigella_flexneri_2a_2457T.faa")
    nucL = fasta.load("seq/ncbi/Shigella_flexneri_2a_2457T.ffn")
    writeSimplifiedAA(aaL,nucL,"seq/ncbi/Shigella_flexneri_2a_2457T.simp.faa")
