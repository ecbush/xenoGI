from Bio import SeqIO
import sys

if __name__ == "__main__":
    outFileName = sys.argv[1]
    outFile = open(outFileName, 'w')
    genbankFileList = sys.argv[2:]
    
    # iterate through list of genbank files
    for fileName in genbankFileList:
        file = open(fileName, 'rU')
        outFile.write(fileName.split("/")[-1][:-5] + "\t")
        # iterate through chromosomes in the genbank file
        for record in SeqIO.parse(file, "genbank"):
            startCodonDict = {}
            #iterate through the genes on the crhomosome
            for feature in record.features:
                # choose only the features that are protein coding genes
                if feature.type == "CDS" and 'protein_id' in feature.qualifiers:
                    start = int(feature.location.start)
                    key = feature.qualifiers['protein_id'][0]
                    startCodonDict[key] = start
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

# keyValL.sort(key=lambda x: x[1])
# geneL = [gene for gene,st in keyValL]
# " ".join(geneL)
