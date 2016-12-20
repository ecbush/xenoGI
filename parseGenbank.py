import sys, glob, os
import genbank

if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))
    
    genbankFileList=glob.glob(params.genbankFilePath)

    # make directory for fasta output
    fastaDir = params.fastaFilePath.split('*')[0]
    os.mkdir(fastaDir)

    # parse
    genbank.parseGenbank(params.geneOrderFN,params.redundProtsFN,params.geneDescriptionsFN,fastaDir,genbankFileList)
