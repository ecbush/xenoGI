import sys
import scores,genomes


if __name__ == "__main__":

    blastFnL=[fn.rstrip() for fn in open(sys.argv[1])]
    protFnL=[fn.rstrip() for fn in open(sys.argv[2])]    

    seqD=genomes.loadProt(protFnL)

    doneSet = set()
    for fn in blastFnL:
        scores.globAlignBlast(fn,seqD,doneSet)


