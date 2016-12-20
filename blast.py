import os, subprocess, glob
from multiprocessing import Pool

def formatDb(dbFileL):
    '''Format fasta files in fastaDir for blast by calling the formatdb
executable.'''
    for dbFileName in dbFileL:
        subprocess.call(['makeblastdb', '-dbtype' ,'prot', '-in', dbFileName],stdout=subprocess.PIPE)
    return

# blastp -matrix BLOSUM62 -gapopen 11 -gapextend 1 -evalue 0.01 -seg yes -outfmt 6 -db fasta/Citrobacter.fa -query fasta/Klebsiella.fa

def makeBlastClineList(dbFileL,fastaFilePath,blastFilePath,blastCLine):
    '''Create a list of lists, where the sublists have the command line
needed to run blastp on a pair of databases.'''

    blastCLineT = tuple(blastCLine.split())

    blastDir = blastFilePath.split("*")[0]
    blastExtension = blastFilePath.split("*")[1]
    
    suffixToRemove = fastaFilePath.split('*')[-1]
    filePathStemToRemove = fastaFilePath.split("*")[0]

    clineL=[]
    for query in dbFileL:
        for db in dbFileL:
            if query != db:
                qstem=query.split(suffixToRemove)[0]
                qstem=qstem.split(filePathStemToRemove)[1]

                dbstem=db.split(suffixToRemove)[0]
                dbstem=dbstem.split(filePathStemToRemove)[1]

                outFN = blastDir + qstem + '-' + dbstem + blastExtension
                L = list(blastCLineT) + ['-query',query,'-db',db,'-out',outFN]
                clineL.append(L)

    return clineL

def subprocessWrapper(cline):
    '''A wraper within which we call subprocess. Becasue we're having
blast write to file, we dump std_out. return std_err.'''
    pipes=subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = pipes.communicate()
    return stderr

def runBlast(fastaFilePath,blastFilePath,blastCLine,numThreads):
    '''Run blast comparing every database against every other in
fastaFilePath. Save to the directory indicated by blastFilePath, using
the blast parameters in blastCLine.'''

    # format the databases
    dbFileL=glob.glob(fastaFilePath)
    formatDb(dbFileL)

    # create blast directory
    blastDir = blastFilePath.split("*")[0]
    os.mkdir(blastDir)
    
    clineL =  makeBlastClineList(dbFileL,fastaFilePath,blastFilePath,blastCLine)

    p=Pool(numThreads)
    stderrL = p.map(subprocessWrapper, clineL)
    
    return
