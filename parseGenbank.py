import sys, glob, os
import parameters,genbank

if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    
    genbankFileList=glob.glob(paramD['genbankFilePath'])

    # if directory for fastas doesn't exist yet, make it
    fastaDir = paramD['fastaFilePath'].split('*')[0]
    if glob.glob(fastaDir)==[]:
        os.mkdir(fastaDir)

    fileNameMapD = parameters.loadFileNameMapD(paramD['fileNameMapFN'])

    # parse
    genbank.parseGenbank(paramD['geneOrderFN'],paramD['redundProtsFN'],paramD['geneInfoFN'],fastaDir,genbankFileList,fileNameMapD)
