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
    

def parseGenbank(geneOrderOutFileName,redundancyOutFileName,geneDescriptionsOutFileName,fastaOutFileDir,genbankFileList):
    '''We pass through each genbank file twice. Once to identify redundant
genes, so we can avoid them. And another time to get the stuff we
want.'''
    
    geneOrderOutFile = open(geneOrderOutFileName, 'w')
    redundFile = open(redundancyOutFileName, 'w')
    geneDescriptionsFile = open(geneDescriptionsOutFileName, 'w')
    
    # iterate through list of genbank files
    for fileName in genbankFileList:

        speciesName = fileName.split("/")[-1][:-5]
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
            #iterate through the genes on the chromosome
            for feature in record.features:
                # choose only the features that are protein coding genes
                if feature.type == "CDS" and 'protein_id' in feature.qualifiers:
                    geneName = speciesName + '-' + feature.qualifiers['protein_id'][0]
                    # verify not in set of genes that appear more than once
                    if geneName in uniqueS:
                        start = int(feature.location.start)
                        aaSeq = feature.qualifiers['translation'][0]

                        # get description
                        if 'product' in feature.qualifiers:
                            descrip = feature.qualifiers['product'][0]
                        else:
                            descrip = ''

                        geneStartSeqL.append((geneName,start,aaSeq,descrip))

            if geneStartSeqL != []: # if its not empty
                # sort by start position
                geneStartSeqL.sort(key=lambda x: x[1])
                geneL=[]
                for geneName,start,aaSeq,descrip in geneStartSeqL:
                    # write to fastaOutFile
                    fastaOutFile.write(">" + geneName + "\n" + aaSeq + "\n")
                    geneDescriptionsFile.write(geneName + "\t" + descrip + "\n")
                    geneL.append(geneName)
                    
                # write this chromosome to geneOrderFile                    
                geneOrderOutFile.write("\t"+" ".join(geneL))
                
        geneOrderOutFile.write("\n")
                    
        inFile.close()
        fastaOutFile.close()

    geneOrderOutFile.close()
    redundFile.close()
    geneDescriptionsFile.close()
