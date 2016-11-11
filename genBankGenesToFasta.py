from Bio import SeqIO
import sys

if __name__ == "__main__":
    outFileDir = sys.argv[1]
    genbankFileList = sys.argv[2:]

    for fileName in genbankFileList:
        file = open(fileName, 'rU')
        outFileName = outFileDir + "/" + fileName.split("/")[-1][:-5] + ".fa"
        outFile = open(outFileName, 'w')
        for record in SeqIO.parse(file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and 'protein_id' in feature.qualifiers:
                    id = feature.qualifiers['protein_id'][0]
                    aaSeq = feature.qualifiers['translation'][0]
                    outFile.write(">" + id + "\n" + aaSeq + "\n")
        outFile.close()
        file.close()
