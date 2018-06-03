import os
from Bio import SeqIO

def getUniqueRedundSets(fileName,speciesName):
    '''Run through genbank file fileName, get set of unique genes (with no
redundancies), and set with redundancies.'''
    f = open(fileName, 'rU')
    geneL=[]
    for record in SeqIO.parse(f, "genbank"):
        # iterate through the genes on the chromosome
        for feature in record.features:
            # choose only the features that are protein coding genes
            if feature.type == "CDS" and 'protein_id' in feature.qualifiers:
                geneName = speciesName + '-' + feature.qualifiers['protein_id'][0]
                geneL.append(geneName)
                
    f.close()

    # now figure out which ones are unique
    uniqueS = set()
    redundS = set()
    for gene in geneL:
        if geneL.count(gene)>1:
            redundS.add(gene)
        else:
            uniqueS.add(gene)
    return uniqueS,redundS
    

def parseGenbank(geneOrderOutFileName,redundancyOutFileName,geneInfoOutFileName,fastaOutFileDir,genbankFileList,fileNameMapD):
    '''We pass through each genbank file twice. Once to identify redundant
genes, so we can avoid them. And another time to get the stuff we
want.'''

    if genbankFileList == []:
        raise ValueError("List of genbank files to parse is empty.")
    
    geneOrderOutFile = open(geneOrderOutFileName, 'w')
    redundFile = open(redundancyOutFileName, 'w')
    geneInfoFile = open(geneInfoOutFileName, 'w')
    
    # iterate through list of genbank files
    for fileName in genbankFileList:

        genbankName = os.path.split(fileName)[-1]

        speciesName = fileNameMapD[genbankName]

        uniqueS,redundS=getUniqueRedundSets(fileName,speciesName)

        # write the redundant ones for this species to our redund file
        for gene in redundS:
            redundFile.write(gene + "\n")
            
        inFile = open(fileName, 'rU')
        fastaOutName = fastaOutFileDir + speciesName + ".fa"
        fastaOutFile = open(fastaOutName, 'w')

        # start next line in geneOrderOutFile
        geneOrderOutFile.write(speciesName)
        
        # iterate through chromosomes in the genbank file
        for record in SeqIO.parse(inFile, "genbank"):
            geneStartSeqL = []
            chrom = record.id # .id, as opposed to .name includes the version id
            #iterate through the genes on the chromosome
            for feature in record.features:
                # choose only the features that are protein coding genes
                if feature.type == "CDS" and 'protein_id' in feature.qualifiers and 'translation' in feature.qualifiers:
                    geneName = speciesName + '-' + feature.qualifiers['protein_id'][0]
                    # verify not in set of genes that appear more than once
                    if geneName in uniqueS:
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        # strand
                        if feature.location.strand == 1:
                            strand = "+"
                        elif feature.location.strand == -1:
                            strand = "-"
                        else:
                            strand = "?"
                        
                        aaSeq = feature.qualifiers['translation'][0]

                        # get common name
                        commonName=''
                        if 'gene' in feature.qualifiers:
                            commonName=feature.qualifiers['gene'][0]

                        locusTag=''
                        if 'locus_tag' in feature.qualifiers:
                            locusTag = feature.qualifiers['locus_tag'][0]
                        
                            
                        # get description
                        descrip=''
                        if 'gene' in feature.qualifiers:
                            descrip += feature.qualifiers['gene'][0]+' - '
                        if 'product' in feature.qualifiers:
                            descrip += feature.qualifiers['product'][0]

                        geneStartSeqL.append((geneName,commonName,locusTag,descrip,chrom,start,end,strand,aaSeq))

            if geneStartSeqL != []: # if its not empty
                # no need to sort, it's in order in the file
                geneL=[]
                for geneName,commonName,locusTag,descrip,chrom,start,end,strand,aaSeq in geneStartSeqL:
                    # write to fastaOutFile
                    fastaOutFile.write(">" + geneName + "\n" + aaSeq + "\n")
                    geneInfoFile.write("\t".join([geneName,commonName,locusTag,descrip,chrom,str(start),str(end),strand]) + "\n")
                    geneL.append(geneName)
                    
                # write this chromosome to geneOrderFile                    
                geneOrderOutFile.write("\t"+" ".join(geneL))
                
        geneOrderOutFile.write("\n")
                    
        inFile.close()
        fastaOutFile.close()

    geneOrderOutFile.close()
    redundFile.close()
    geneInfoFile.close()
