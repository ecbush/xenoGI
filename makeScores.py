import sys
import scores,genomes


if __name__ == "__main__":

    blastFnL=[fn.rstrip() for fn in open(sys.argv[1])]
    protFnL=[fn.rstrip() for fn in open(sys.argv[2])]    

    paramFN=sys.argv[3]
    params = __import__(paramFN.replace('.py', ''))

    
    seqD=genomes.loadProt(protFnL)

    scores.createSimScores(blastFnL,seqD,params.numThreads,params.scoresFN)
