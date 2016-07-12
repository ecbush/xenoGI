import sys,fasta,parasail

## Funcs

def loadProt(protFnL):
    '''Given a list of file names of simp.faa type fasta files, load the
sequences and store in a dictionary keyed by protein name.
    '''
    seqD={}
    for fn in protFnL:
        for header,seq in fasta.load(fn):
            gn = header.split()[0][1:]
            seqD[gn]=seq
    return seqD

def globAlignBlast(fn,seqD,doneSet):
    '''Given a file name with blast output, go through each hit and run
needleman wunch on the sequences. Print gene names, score and blast
percID.
    '''
    f=open(fn,'r')
    while True:
        s = f.readline()
        if s=='':
            break
        L = s.split('\t')
        if len(L) != 12: # we only want lines with 12 columns
            continue

        g1 = L[0]
        g2 = L[1]
        percID = L[2]

        if (g1,g2) in doneSet: continue

        r=parasail.nw_scan(seqD[g1],seqD[g2], 11, 1, parasail.blosum62)

        print(g1,g2,r.score,percID,sep='\t')

        doneSet.add((g1,g2))
        doneSet.add((g2,g1))

        # note. parasail stats is currently messed up, giving
        # wrong length. We'll just use score here. In future, when
        # that package allows you to get the alignment, we can
        # work on getting true distances.

        # also, clearly we need to parallelize this in future.

    f.close()
    

    
## Main

if __name__ == "__main__":

    blastFnL=[fn.rstrip() for fn in open(sys.argv[1])]
    protFnL=[fn.rstrip() for fn in open(sys.argv[2])]    

    seqD=loadProt(protFnL)

    doneSet = set()
    for fn in blastFnL:
        globAlignBlast(fn,seqD,doneSet)
