import sys
import blast

if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    blast.runBlast(params.fastaFilePath,params.blastFilePath,params.blastCLine,params.numThreads)
