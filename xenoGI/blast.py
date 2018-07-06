import sys, os, subprocess, glob
from . import trees
from multiprocessing import Pool
from Bio import Phylo

def getDbFileL(fastaFilePath,treeFN):
    '''Obtains and returns a list of all fasta files that should be run
through blast. If there is a valid tree file present, then this list
should contain only those fastas that correspond to species in the
tree. If there is not a valid tree file, then the list will contain
all fasta files present in the fasta directory. The reason we don't
want to assume there is a tree file is that sometimes users will need
to run blast before making that.
    '''
    # get list of all fastas in directory
    dbFileL=glob.glob(fastaFilePath)

    # load tree if present
    try:
        newDbFileL=[]
        tree,strainStr2NumD,strainNum2StrD = trees.readTree(treeFN)
        leavesL=[strainNum2StrD[strainNum] for strainNum in trees.leafList(tree)]

        for dbFile in dbFileL:
            dbStem = os.path.split(dbFile)[-1] # strain + extension
            dbStem = os.path.splitext(dbStem)[0] # strain only
            if dbStem in leavesL:
                # this file name is in the tree, include it
                newDbFileL.append(dbFile)

        dbFileL = newDbFileL
    except ( TypeError, FileNotFoundError, Phylo.NewickIO.NewickError ):
        print("Since a properly formated tree file was not provided, we can't use that to determine which fasta files to blast. We will therefore blast all fasta files in fastaFilePath.",file=sys.stderr)

    return dbFileL

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

def runBlast(fastaFilePath,blastFilePath,blastExecutDirPath,blastCLine,numThreads,treeFN):
    '''Run blast comparing every database against every other in
fastaFilePath. Save to the directory indicated by blastFilePath, using
the blast parameters in blastCLine.'''

    # format the databases
    dbFileL=getDbFileL(fastaFilePath,treeFN)

    formatDb(dbFileL,blastExecutDirPath)
    
    # if directory for blast doesn't exist yet, make it
    blastDir = os.path.split(blastFilePath)[0]
    if glob.glob(blastDir)==[]:
        os.mkdir(blastDir)

    clineL =  makeBlastClineList(dbFileL,fastaFilePath,blastFilePath,blastExecutDirPath,blastCLine)

    p=Pool(numThreads)
    stderrL = p.map(subprocessWrapper, clineL)
    
    return
