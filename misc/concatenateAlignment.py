import sys,os
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import fasta

# We assume each sequence is on a single line

def concatenateAlign(inputFileList, outputFile):
    '''Reads a file of aligned genes, and creates a super alignment''' 

    concatAlignmentDict = {} #Dict will be keyed on strain name, and values will be a list.
    
    for fn in inputFileList:

        for hd,seq in fasta.load(fn):
            strainName = hd.split(' ')[0]

            if strainName in concatAlignmentDict:
                concatAlignmentDict[strainName].append(seq)
            else:
                concatAlignmentDict[strainName] = [seq]


    # The dict should now contain everything. 

    fOut = open(outputFile, "w")

    for species in concatAlignmentDict:

        fOut.writelines([species+"\n", "".join(concatAlignmentDict[species])+"\n"])

    fOut.close()


if __name__ == "__main__":

        outputFile = sys.argv[1]
        inputFileList = sys.argv[2:]
        
        concatenateAlign(inputFileList, outputFile)
