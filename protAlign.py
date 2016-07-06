# EB
import sys, fasta, os

#Usage: python alignProt.py full/orthologs.out full/seqLists.txt full/prot.afa

## funcs

def loadOrthos(orthoFN):
    '''Reads the orthologs file. The first lines of the file contains
species names, in the form: number species name. Gets all the names in
a list. Then reads the remaining lines to get ortholog sets. Puts each
set in a tuple, collecting the tuples in a list. Returns these two
lists.

    '''
    f= open(orthoFN,'r')

    strainL=[]
    orthoL=[]
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.rstrip().split()
        if s[0].isdigit(): # is a number-species line
            strainL.append(L[1])
        else:
            # is an ortho list line
            orthoL.append(tuple(L))
            
            
    f.close()
    return strainL,orthoL

def loadSeqs(strainSeqFileListFN):
    '''Load list of protein sequence files, then load each of these and
put in a dict (sequence dict, keyed by gene name). Keep these dicts in
a dict keyed by strain name

    '''
    seqFileL=[x.rstrip() for x in open(strainSeqFileListFN,"r")]

    seqD={} # will be keyed by strain name
    for seqFile in seqFileL:
        D={} # will be keyed by gene
        seqL=fasta.load(seqFile)
        for hd,sq in seqL:
            gnNm=hd[1:].split()[0] # lop off >, then split on white space
            D[gnNm]=sq
        strainName = seqFile.split(".simp.faa")[0].split("/")[2]
        seqD[strainName]=D
    return seqD

def writeSeqBlock(f,strainL,orthos,seqD):
    '''writes a multifasta block. seqs are specified in
orthos, and obtained from seqD. the order of species in orthos is
given by strainL. And we must use the strain name to get the right
subdict from seqD.'''
    for i in range(len(orthos)):
        f.write(">"+orthos[i]+"\n")
        f.write(seqD[strainL[i]][orthos[i]]+"\n")
    f.write("\n") # blank line between blocks
        
if __name__ == "__main__":

   orthoFN = sys.argv[1]
   strainSeqFileListFN = sys.argv[2]
   alignmentFN = sys.argv[3]

   strainL,orthoL=loadOrthos(orthoFN)
   seqD=loadSeqs(strainSeqFileListFN)

   # write these to alignmentFN
   f=open(alignmentFN,"w")

   intempAlignFN="/tmp/tempAlign.fa"
   outtempAlignFN="/tmp/tempAlign.afa"
   
   for orthos in orthoL:
       tempf = open(intempAlignFN,"w")
       writeSeqBlock(tempf,strainL,orthos,seqD)
       tempf.close()

       # align the temp file
       os.system("muscle -in "+ intempAlignFN + " -out " + outtempAlignFN)

       # load aligned file into a sequence dict
       alSeqL=fasta.load(outtempAlignFN)
       alSeqD={}
       for hd,sq in alSeqL:
           alSeqD[hd[1:]]=sq

       # write aligned fasta block into main output file, in the
       # original order (muscle messes up this order).
       for gene in orthos:
           f.write(">"+gene+"\n")
           f.write(alSeqD[gene]+"\n")
       f.write("\n")
       
   f.close()
