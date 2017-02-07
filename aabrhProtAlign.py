import sys,glob,os,random
import genomes,fasta,parameters

## funcs

def loadOrthos(aabrhFN):
    '''Reads the all around best reciprocal hits orthologs file. One set
per line.'''
    f= open(aabrhFN,'r')

    orthoL=[]
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.rstrip().split('\t')
        orthoL.append(tuple(L))
            
    f.close()
    return orthoL


def writeSeqBlock(f,orthos,seqD):
    '''writes a multifasta block. seqs are specified in
orthos, and obtained from seqD. the order of species in orthos is
given by strainL. And we must use the strain name to get the right
subdict from seqD.'''
    for i in range(len(orthos)):
        f.write(">"+orthos[i]+"\n")
        f.write(seqD[orthos[i]]+"\n")
    f.write("\n") # blank line between blocks
     

if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)


    aabrhAlignmentFN = sys.argv[2]
    
    aabrhL = loadOrthos(paramD['aabrhFN'])
    protFnL=glob.glob(paramD['fastaFilePath'])
    seqD=genomes.loadProt(protFnL)
    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])
    
    # write these to alignmentFN
    f=open(aabrhAlignmentFN,"w")

    randStr = str(random.randrange(1e5))
    intempAlignFN="/tmp/tempAlign"+randStr+".fa"
    outtempAlignFN="/tmp/tempAlign"+randStr+".afa"
   
    for orthos in aabrhL:
        tempf = open(intempAlignFN,"w")
        writeSeqBlock(tempf,orthos,seqD)
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
            commonName,locusTag,descrip,chrom,start,end,strand=geneInfoD[gene]
            f.write(">"+gene+" "+locusTag+"\n")
            f.write(alSeqD[gene]+"\n")
        f.write("\n")

    # delete the temp files
    os.system("rm "+intempAlignFN)
    os.system("rm "+outtempAlignFN)
    
    f.close()

