import sys,os,glob
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import parameters,trees,genomes,scores


if __name__ == "__main__":

    paramFN=sys.argv[1]
    dnaOrProt = sys.argv[2]
    orthosFN = sys.argv[3]
    alignDirName = sys.argv[4]

    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)
    genesO = genomes.genes(paramD['geneInfoFN'])

    # load seqs
    protSeqD=genomes.loadSeq(paramD, '_prot.fa')
    if dnaOrProt == "dna":
        dnaSeqD = genomes.loadSeq(paramD, '_dna.fa')
    else:
        dnaSeqD = {}

    # load ortho groups
    aabrhHardCoreL = scores.loadOrthos(orthosFN)

    # create alignDir if it doesn't already exist
    if glob.glob(alignDirName)==[]:
        os.mkdir(alignDirName)
    
    # align
    inTempProtFN=os.path.join(alignDirName,"tempProt.fa")
    orthoGroupNum = 0
    for orthoT in aabrhHardCoreL:
        orthoGroupNumStr = str(orthoGroupNum).zfill(6) # pad with 0's so ls will display in right order
        outAlignFN=os.path.join(alignDirName,"align"+orthoGroupNumStr+".afa")
        trees.alignOneOrthoT(orthoT,True,paramD['musclePath'],inTempProtFN,outAlignFN,protSeqD,dnaSeqD,genesO)
        orthoGroupNum+=1

    os.remove(inTempProtFN)
