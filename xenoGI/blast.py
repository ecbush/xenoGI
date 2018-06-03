import os, subprocess, glob
from multiprocessing import Pool

def formatDb(dbFileL,blastExecutDirPath):
    '''Format fasta files in fastaDir for blast by calling the formatdb
executable.'''
    for dbFileName in dbFileL:
        makeblastdbExecutable = os.path.join(blastExecutDirPath,'makeblastdb')
        subprocess.call([makeblastdbExecutable, '-dbtype' ,'prot', '-in', dbFileName],stdout=subprocess.PIPE)
    return

def makeBlastClineList(dbFileL,fastaFilePath,blastFilePath,blastExecutDirPath,blastCLine):
    '''Create a list of lists, where the sublists have the command line
needed to run blastp on a pair of databases.'''

    # get a tuple of the blastp command line args for use below
    blastCLineL = blastCLine.split()
    blastCLineL[0] = os.path.join(blastExecutDirPath,blastCLineL[0])
    blastCLineT = tuple(blastCLineL)

    # get blast file dir and extension
    splitT = os.path.split(blastFilePath)
    blastDir,rest = splitT
    blastExtension = rest.split("*")[-1]

    clineL=[]
    for query in dbFileL:
        for db in dbFileL:
            if query != db:

                qstem = os.path.split(query)[-1]
                qstem = os.path.splitext(qstem)[0]
                
                dbstem = os.path.split(db)[-1]
                dbstem = os.path.splitext(dbstem)[0]

                outFN = os.path.join( blastDir, qstem + '-' + dbstem + blastExtension )
                L = list(blastCLineT) + ['-query',query,'-db',db,'-out',outFN]
                clineL.append(L)

    return clineL

def subprocessWrapper(cline):
    '''A wraper within which we call subprocess. Becasue we're having
blast write to file, we dump std_out. return std_err.'''
    pipes=subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = pipes.communicate()
    return stderr

def runBlast(fastaFilePath,blastFilePath,blastExecutDirPath,blastCLine,numThreads):
    '''Run blast comparing every database against every other in
fastaFilePath. Save to the directory indicated by blastFilePath, using
the blast parameters in blastCLine.'''

    # format the databases
    dbFileL=glob.glob(fastaFilePath)
    formatDb(dbFileL,blastExecutDirPath)
    
    # if directory for blast doesn't exist yet, make it
    blastDir = os.path.split(blastFilePath)[0]
    if glob.glob(blastDir)==[]:
        os.mkdir(blastDir)

    clineL =  makeBlastClineList(dbFileL,fastaFilePath,blastFilePath,blastExecutDirPath,blastCLine)

    p=Pool(numThreads)
    stderrL = p.map(subprocessWrapper, clineL)
    
    return
