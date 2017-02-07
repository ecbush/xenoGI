import sys
import blast,parameters

if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    blast.runBlast(paramD['fastaFilePath'],paramD['blastFilePath'],paramD['blastCLine'],paramD['numThreads'])
