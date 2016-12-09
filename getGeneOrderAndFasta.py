from Bio import SeqIO
import sys

if __name__ == "__main__":
    outFileName = sys.argv[1] #name for gene order file
    redundancyFileName = sys.argv[2] #name for file of redundant proteins
    fastaFileDir = sys.argv[3] #relative directory to send fasta files
    outFile = open(outFileName, 'w')
    redundFile = open(redundancyFileName, 'w')
    genbankFileList = sys.argv[4:]

    #initialize dictionary to keep track of proteins we've seen already
    seenDict = {}
    
    # iterate through list of genbank files
    for fileName in genbankFileList:
        file = open(fileName, 'rU')
        speciesName = fileName.split("/")[-1][:-5]
        fastaOutName = fastaFileDir + speciesName + ".fa"
        fastaFile = open(fastaOutName, 'w')
        outFile.write(speciesName + "\t")
        redundFile.write("\n" + speciesName + ":" + "\n")
        # iterate through chromosomes in the genbank file
        for record in SeqIO.parse(file, "genbank"):
            startCodonDict = {}
            #iterate through the genes on the crhomosome
            for feature in record.features:
                # choose only the features that are protein coding genes
                if feature.type == "CDS" and 'protein_id' in feature.qualifiers:
                    key = speciesName + feature.qualifiers['protein_id'][0]
                    #check if the protein has been seen before
                    if key not in seenDict:
                        start = int(feature.location.start)
                        #key = feature.qualifiers['protein_id'][0]
                        startCodonDict[key] = start
                        seenDict[key] = True
                        aaSeq = feature.qualifiers['translation'][0]
                        fastaFile.write(">" + key + "\n" + aaSeq + "\n")
                    #if it has, write it to the redundancy file and remove it from the dict
                    else:
                        if key in startCodonDict:
                            del startCodonDict[key]
                            redundFile.write(key + "\n")
            # make a list of entries in the dict and sort them by start position
            keyValL = list(startCodonDict.items())
            # make a string with all the gene names in order of start position
            geneL = [gene for gene,str in keyValL]
            # write the string to the out file
            outFile.write(" ".join(geneL))
            outFile.write("\t")
        file.close()
        outFile.write("\n")
    outFile.close()
