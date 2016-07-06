import sys

# general idea. take a list of broad summary files and a second list of ncbi. use to load and get gene order.


#LOCUS   SYMBOL  SYNOYM  LENGTH  START   STOP    STRAND  NAME    CHROMOSOME      GENOME ONTOLOGY ENZYME CODE     KEGG    PATHWAY REACTION        



# Broad

def loadBroad(fn):
    '''Read in the broad genome summary file, and return a dictionary
keyed by chromosome number with a sorted list of (start,LOCUS) tuples
as value.'''
    D={}
    f=open(fn,'r')
    s=f.readline() # skip header line
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.split('\t')
        locus=L[0]
        start=int(L[4])
        chrom=L[8]

        if chrom in D:
            D[chrom].append((start,locus))
        else:
            D[chrom]=[(start,locus)]

    f.close()
    
    for key in D:
        D[key].sort()

    # get strain name from fn
    strainName=fn.split('genome_summary_per_gene.txt')[0].split('/')[-1]
        
    return strainName,D

def mkGeneOrderStr(strainName,D):
    '''Return a string with strain name. then contig units separated by |. and genes within a contig separated by ' '.'''
    contigL=[]
    for key in D:
        contigL.append(' '.join((geneName for st,geneName in D[key])))

    return strainName + '\t' + '\t'.join(contigL)

def loadNCBI(fn):
    '''Given NCBI data in our .simp.faa format, get the gene names and start locations. Should replace this with something better (something that uses .gbff files) in future.'''


    D={}
    D['chr']=[] # always this name for chr here, since we know the ones we're working with for now only have one chr. will need to change in future, when we use more ncbi stuff.

    f=open(fn,'r')
    while True:
        h=f.readline()
        if h=='':break
        if h[0]=='>':

            j=h.split()[1]
            j=j.split('-')[0]
            if j[0]=='c':
                j=j[1:]
            j=int(j)
            name=h.split()[0][1:]
            D['chr'].append((j,name))
    f.close()

    D['chr'].sort()
    
    # get strain name from fn
    strainName=fn.split('.simp.faa')[0].split('/')[-1]
    
    return strainName,D


## Main

if __name__ == "__main__":

    broadFnL=[fn.rstrip() for fn in open(sys.argv[1])]
    ncbiFnL=[fn.rstrip() for fn in open(sys.argv[2])]

    for fn in broadFnL:
        print(mkGeneOrderStr(*loadBroad(fn)))

    for fn in ncbiFnL:
        print(mkGeneOrderStr(*loadNCBI(fn)))
