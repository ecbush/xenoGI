import sys
import parameters,blast

if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    blast.runBlast(paramD['fastaFilePath'],paramD['blastFilePath'],paramD['blastExecutDirPath'],paramD['blastCLine'],paramD['numThreads'])
