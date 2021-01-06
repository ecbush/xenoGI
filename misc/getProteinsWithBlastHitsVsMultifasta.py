import sys,os
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import parameters,blast,genomes,xenoGI

# Identifies xenoGI proteins with significant similarity to a provided
# protein multifasta, printing their xenoGI gene names to
# standard out, one per line. Expects to be run in a xenoGI working
# directory.

if __name__ == "__main__":

    paramFN = sys.argv[1]
    strainInfoFN = sys.argv[2] # include as argument so can pass in scaffold only set in xlMode case
    proteinMultiFastaPath = sys.argv[3]

    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)
    genesO = genomes.genes(paramD['geneInfoFN'])
    genesO.initializeGeneNumToNameD(paramD['geneInfoFN'])
    strainNamesT = xenoGI.readStrainInfoFN(strainInfoFN)
    
    evalueThresh = paramD['evalueThresh']
    alignCoverThresh = paramD['alignCoverThresh']
    percIdentThresh =  paramD['percIdentThresh']
    blastFilePath = paramD['blastFilePath']
    blastFileJoinStr = paramD['blastFileJoinStr']

    # Blast them vs. all strains
    fastaDir = paramD['fastaFilePath'].split('*')[0]
    allStrainsFileNamesL = []
    for strain in strainNamesT:
        allStrainsFileNamesL.append(os.path.join(fastaDir,strain+"_prot.fa"))

    blast.runBlast([proteinMultiFastaPath],allStrainsFileNamesL,paramD)
    
    # parse blast to identify ifams with similarity to model seqs
    blastDir = os.path.split(blastFilePath)[0]

    proteinMultiFastaFN = os.path.split(proteinMultiFastaPath)[-1]
    # remove _prot.fa
    proteinMultiFastaFN = proteinMultiFastaFN.split("_prot.fa")[0]

    # if additional extension
    proteinMultiFastaStem = os.path.splitext(proteinMultiFastaFN)[0]
    
    outS = set()
    for strainName in strainNamesT:
        # remove dir (if any) from proteinMultiFastaPath
        fileStr = proteinMultiFastaStem+blastFileJoinStr+strainName+'.out'
        fn = os.path.join(blastDir,fileStr)
        for g1,g2,evalue,alCov,pident,score in blast.parseBlastFile(fn,evalueThresh,alignCoverThresh,percIdentThresh):    
            # g2 will always be a xenoGI gene
            outS.add((g2,genesO.numToName(g2)))
            
    # sort by gene number and print
    for geneNum,geneStr in sorted(outS):
        print(geneStr)

