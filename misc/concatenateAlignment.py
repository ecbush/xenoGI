import sys

# We assume each sequence is on a single line

def concatenateAlign(alignedBlockFile, outputFile):
    '''Reads a file of aligned genes, and creates a super alignment''' 

    f=open(alignedBlockFile,'r')

    concatAlignmentDict = {} #Dict will be keyed on strain name, and values will be a list.

    while True:

        s = f.readline()
        if s == '':
            break
        
        elif s[0] == ">":
            strainName = s.split(' ')[0] #Takes a line, splits on the empty space, then takes the first part of the list and keeps the part with the strain name.
            strainName = strainName.split("-")[0] #Save the strain name, get rid of the part after the "-", which is protein specific
    
        elif s[0] == "\n":
            pass 
    
        else: #This should be the sequence, so add it to the dict
            sequence = s.rstrip() #rstrip removes trailing whitespace

            if strainName in concatAlignmentDict:
                concatAlignmentDict[strainName].append(sequence)
            else:
                concatAlignmentDict[strainName] = [sequence]

    f.close()

    # The dict should now contain everything. 

    fOut = open(outputFile, "w")

    for species in concatAlignmentDict:

        fOut.writelines([species+"\n", "".join(concatAlignmentDict[species])+"\n"])

    fOut.close()


if __name__ == "__main__":

        alignedBlockFile = sys.argv[1]
        outputFile = sys.argv[2]

        concatenateAlign(alignedBlockFile, outputFile)
