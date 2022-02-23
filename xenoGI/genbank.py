import os,sys
from Bio import SeqIO

def parseGenbank(paramD,fastaOutFileDir,genbankFileList,fileNameMapD):
    '''Parse all the genbank (gbff) files in genbankFileList.'''
    
    if genbankFileList == []:
        raise ValueError("List of genbank files to parse is empty.")

    dnaBasedGeneTrees = paramD['dnaBasedGeneTrees']
    
    geneInfoFile = open(paramD['geneInfoFN'], 'w')
    geneOrderOutFile = open(paramD['geneOrderFN'], 'w')
    
    problemGenbankFileL = []
    
    # iterate through list of genbank files
    geneNum =  0 # xenoGI internal gene numbering
    for fileName in genbankFileList:

        geneNum,problemGenbankFileL = parseGenbankSingleFile(geneNum,fileName,dnaBasedGeneTrees,fileNameMapD,geneInfoFile,geneOrderOutFile,fastaOutFileDir,problemGenbankFileL)
        
    geneInfoFile.close()
    geneOrderOutFile.close()

    # If there are any files in problemGenbankFileL, throw error
    if problemGenbankFileL != []:

        with open(paramD['problemGenbankFN'],'w') as problemGenbankF:
            for probT in problemGenbankFileL:
                problemGenbankF.write("\t".join(probT)+"\n")
    
        raise ValueError('Some genbank files have problems with their annotations. They are listed in ' + paramD['problemGenbankFN'] + '. Please remove and run again.\n')
    
def parseGenbankSingleFile(geneNum,fileName,dnaBasedGeneTrees,fileNameMapD,geneInfoFile,geneOrderOutFile,fastaOutFileDir,problemGenbankFileL):
    '''Parse a single genbank file.
    '''
    genbankName = os.path.split(fileName)[-1]
    speciesName = fileNameMapD[genbankName]

    # verify that there are protein annotations
    if not verifyProteinAnnotations(fileName):
        problemGenbankFileL.append((fileName,'lacks protein annotations'))
        return geneNum,problemGenbankFileL 

    # if using dna for gene trees, verify that dna annotations 3x
    # longer than protein
    if dnaBasedGeneTrees and not verifyDnaAnnotations(fileName):
        problemGenbankFileL.append((fileName,'length of dna and protein annotations do not correspond properly'))
        return geneNum,problemGenbankFileL 
    
    # start a block for this species in geneInfoFile
    geneInfoFile.write("# "+speciesName+"\n")
    
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
            if feature.type == "CDS" and 'translation' in feature.qualifiers:

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

                # common name
                commonName=''
                if 'gene' in feature.qualifiers:
                    commonName=feature.qualifiers['gene'][0]

                # locus tag
                locusTag=''
                if 'locus_tag' in feature.qualifiers:
                    locusTag = feature.qualifiers['locus_tag'][0]

                # protein ID
                proteinId=''
                if 'protein_id' in feature.qualifiers:
                    proteinId = feature.qualifiers['protein_id'][0]

                # get description
                descrip=''
                if 'gene' in feature.qualifiers:
                    descrip += feature.qualifiers['gene'][0]+' - '
                if 'product' in feature.qualifiers:
                    descrip += feature.qualifiers['product'][0]

                # get xenoGI gene name
                if locusTag != '':
                    geneName = str(geneNum) + "_" + speciesName + '-' + locusTag
                else:
                    # if locus tag missing, base name on start and end coordinates
                    geneName = str(geneNum) + "_" + speciesName + '-' + str(start)+"_"+str(end)
                    
                # write to fastaOutFile
                protFastaOutFile.write(">" + geneName + "\n" + aaSeq + "\n")
                if dnaBasedGeneTrees:
                    dnaFastaOutFile.write(">" + geneName + "\n" + dnaSeq + "\n")


                # write to gene information file
                geneInfoFile.write("\t".join([str(geneNum),geneName,commonName,locusTag,proteinId,descrip,chrom,str(start),str(end),strand]) + "\n")
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

def verifyProteinAnnotations(fileName):
    '''Run through genbank file fileName verifying that it contains
protein annotations. Return True if it does, False if not.'''
    with open(fileName, 'rU') as f:
        for record in SeqIO.parse(f, "genbank"):
            # iterate through the genes on the chromosome
            for feature in record.features:
                # choose only the features that are protein coding genes
                if feature.type == "CDS" and 'translation' in feature.qualifiers:
                    # at least one annotation
                    return True
    # made it all the way through with no annotations
    return False
    
def verifyDnaAnnotations(fileName):
    '''Run through genbank file fileName verifying that dna gene
annotations are 3x longer than protein annotations. Return True if so
in all cases, False otherwise.
    '''
    with open(fileName, 'rU') as f:
        for record in SeqIO.parse(f, "genbank"):
            # iterate through the genes on the chromosome
            for feature in record.features:
                # choose only the features that are protein coding genes
                if feature.type == "CDS" and 'protein_id' in feature.qualifiers and 'translation' in feature.qualifiers:
                    aaSeq = feature.qualifiers['translation'][0]
                    dnaSeq = str(feature.extract(record.seq))
                    if (len(aaSeq) + 1) * 3 != len(dnaSeq):
                        return False
    # made it all the way through so annotations ok
    return True
   
