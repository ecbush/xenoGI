import sys
import genomes,scores

if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    ## load data structures we'll use below
    geneNames = genomes.geneNames(params.geneOrderFN)
    geneOrderD = genomes.createGeneOrderD(params.geneOrderFN,geneNames)
    strainNamesL = sorted(geneOrderD.keys())


    aabrhL = scores.createAabrhL(params.blastFilePath,strainNamesL,params.evalueThresh,params.aabrhFN)
