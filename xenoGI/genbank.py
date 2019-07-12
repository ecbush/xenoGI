import os,sys
from Bio import SeqIO

def parseGenbank(paramD,fastaOutFileDir,genbankFileList,fileNameMapD):
    '''Parse all the genbank (gbff) files in genbankFileList.'''
    
    if genbankFileList == []:
        raise ValueError("List of genbank files to parse is empty.")

    dnaBasedGeneTrees = paramD['dnaBasedGeneTrees']
    
    geneInfoFile = open(paramD['geneInfoFN'], 'w')
    redundFile = open(paramD['redundProtsFN'], 'w')
    geneOrderOutFile = open(paramD['geneOrderFN'], 'w')
    
    problemGenbankFileL = []
    
    # iterate through list of genbank files
    geneNum =  0 # xenoGI internal gene numbering
    for fileName in genbankFileList:

        geneNum,problemGenbankFileL = parseGenbankSingleFile(geneNum,fileName,dnaBasedGeneTrees,fileNameMapD,geneInfoFile,redundFile,geneOrderOutFile,fastaOutFileDir,problemGenbankFileL)
        
    geneInfoFile.close()
    redundFile.close()
    geneOrderOutFile.close()

    # If there are any files in problemGenbankFileL, throw error
    if problemGenbankFileL != []:
    
        raise ValueError("The following genbank files lack protein annotations:\n" + "\n".join(problemGenbankFileL))
    
def parseGenbankSingleFile(geneNum,fileName,dnaBasedGeneTrees,fileNameMapD,geneInfoFile,redundFile,geneOrderOutFile,fastaOutFileDir,problemGenbankFileL):
    '''Parse a single genbank file. We pass through twice. Once to
identify redundant genes, so we can avoid them. And another time to
get the stuff we want.

    '''
    genbankName = os.path.split(fileName)[-1]
    speciesName = fileNameMapD[genbankName]

    # start a block for this species in geneInfoFile
    geneInfoFile.write("# "+speciesName+"\n")

    uniqueS,redundS=getUniqueRedundSets(fileName,speciesName)

    # Check if we've been passed a .gbff file with no protein
    # annotations, and if so tell user
    if len(uniqueS) == 0 and len(redundS) == 0:
        problemGenbankFileL.append(fileName)
        return geneNum,problemGenbankFileL 

    # write the redundant ones for this species to our redund file
    for gene in redundS:
        redundFile.write(gene + "\n")

    inFile = open(fileName, 'rU')
    protFastaOutName = fastaOutFileDir + speciesName + "_prot.fa"
    protFastaOutFile = open(protFastaOutName, 'w')

    if dnaBasedGeneTrees:
        dnaFastaOutName = fastaOutFileDir + speciesName + '_dna.fa'
        dnaFastaOutFile = open(dnaFastaOutName, 'w')
        
    # start next line in geneOrderOutFile
    geneOrderOutFile.write(speciesName)

    # iterate through chromosomes in the genbank file
    for record in SeqIO.parse(inFile, "genbank"):
        chrom = record.id # .id, as opposed to .name includes the version id
        # iterate through the genes on the chromosome
        genesOnChromL = []
        for feature in record.features:
            # choose only the features that are protein coding genes
            if feature.type == "CDS" and 'protein_id' in feature.qualifiers and 'translation' in feature.qualifiers:
                speciesProtName = speciesName + '-' + feature.qualifiers['protein_id'][0]
                # verify not in set of genes that appear more than once
                if speciesProtName in uniqueS:
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
                    if dnaBasedGeneTrees:
                        dnaSeq = str(feature.extract(record.seq))

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

                    # write to fastaOutFile
                    geneName = str(geneNum) + "_" + speciesProtName
                    

                    protFastaOutFile.write(">" + geneName + "\n" + aaSeq + "\n")
                    if dnaBasedGeneTrees:
                        dnaFastaOutFile.write(">" + geneName + "\n" + dnaSeq + "\n")

                    
                    # write to gene information file
                    geneInfoFile.write("\t".join([str(geneNum),geneName,commonName,locusTag,descrip,chrom,str(start),str(end),strand]) + "\n")
                    genesOnChromL.append(str(geneNum))

                    geneNum += 1 # increment every time we write a line to geneInfoFile
                    

        if genesOnChromL != []:
            # if not empty, write this chromosome to geneOrderFile                    
            geneOrderOutFile.write("\t"+" ".join(genesOnChromL))

    geneOrderOutFile.write("\n") # add newline to geneOrder file, as we're done with this strain.

    inFile.close()
    protFastaOutFile.close()
    if dnaBasedGeneTrees:
        dnaFastaOutFile.close()

    return geneNum,problemGenbankFileL 

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
                speciesProtName = speciesName + '-' + feature.qualifiers['protein_id'][0]
                geneL.append(speciesProtName)
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
