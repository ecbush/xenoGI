import sys,fasta,parasail

def simScore(s1,s2):
    '''Calculate score between a pair of protein sequences, based on a
global alignment. We scale the alignment score to be between 0 and 1,
based on the max and min possible scores for these sequences..'''

    opn = 12
    ext = 1
    matr = parasail.blosum62
    # note, since parasail doesn't charge extend on the first base its
    # like emboss, and not like ncbi blast. NCBI blast uses existance 11,
    # extend 1 for blosum 62, thus we should use open 12 ext 1 here.

    
    r_s1s2 = parasail.nw_scan(s1,s2, opn, ext, matr)

    if len(s1) < len(s2):
        r_self = parasail.nw_scan(s1,s1, opn, ext, matr)
    else:
        r_self = parasail.nw_scan(s2,s2, opn, ext, matr)

    sc = r_s1s2.score
    mx = r_self.score # max possible is shorter seq against itself.
    
    # lowest possible score, if we have gaps opposite all residues and
    # two opens. Note parasail does not count gap extend for the
    # residue where a gap is opened, hence the -2 in the extend
    # formula below.
    mn = - ( 2 * opn + ( (len(s1)+len(s2) -2 ) * ext ) )

    scaled = (sc - mn) / (mx - mn)
    
    return scaled
    
def globAlignBlast(fn,seqD,doneSet):
    '''Given a file name with blast output, go through each hit and run
needleman wunch on the sequences. Print gene names, simScore and blast
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

        scaled = simScore(seqD[g1],seqD[g2])

        print(g1,g2,format(scaled,".6f"),sep='\t')

        doneSet.add((g1,g2))
        doneSet.add((g2,g1))

        # note. parasail stats is currently messed up, giving
        # wrong length. We'll just use score here. In future, when
        # that package allows you to get the alignment, we can
        # work on getting something based on non-gap sites only.

        # also, clearly we need to parallelize this in future.

    f.close()
    
