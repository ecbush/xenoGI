import os,sys
from Bio import SeqIO

fileName = sys.argv[1]

# Usage: python3 ~/xenoGI/misc/examineGenbankDnaVsProtLen.py ncbi/someGenome.gbff
# prints out those cases where dna/3 is the same as aa length
# note that normally dna/3 should equal aa len + 1 (accounting for the stops being missing in the aa seqs).

if __name__ == "__main__":

    with open(fileName, 'rU') as f:
        for record in SeqIO.parse(f, "genbank"):
            # iterate through the genes on the chromosome
            for feature in record.features:
                # choose only the features that are protein coding genes
                if feature.type == "CDS" and 'translation' in feature.qualifiers:
                    aaSeq = feature.qualifiers['translation'][0]
                    dnaSeq = str(feature.extract(record.seq))
                    if (len(aaSeq) + 1) * 3 != len(dnaSeq):
                        print("lengths","aa",len(aaSeq),"dna/3",len(dnaSeq)/3)
                        print(aaSeq)
                        print(dnaSeq)
                        print("-----")
