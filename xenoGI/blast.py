import sys, os, subprocess, glob
from . import trees
from multiprocessing import Pool
from Bio import Phylo

def runBlast(dbFileL_1,dbFileL_2,paramD):
    '''Run blast comparing every database in dbFileL_1 against every
database in dbFileL_2.
    '''

    # if directory for blast doesn't exist yet, make it
    blastDir = os.path.split(paramD['blastFilePath'])[0]
    if glob.glob(blastDir)==[]:
        os.mkdir(blastDir)

    clineL =  makeBlastClineList(dbFileL_1,dbFileL_2,paramD)

    if clineL != []:
        # we need to run some
        
        # format the databases
        uniqueDbL=list(set(dbFileL_1+dbFileL_2))
        formatDb(uniqueDbL,paramD['blastExecutDirPath'])

        with Pool(processes=paramD['numProcesses']) as p:
            for stderr in p.imap_unordered(subprocessWrapper, clineL):
                pass # ignore stderr

            
def getDbFileL(fastaFilePath,strainNamesT):
    '''Obtains and returns a list of all fasta files that should be run
through blast. Only keeps those that are in stainNamesT.
    '''
    # get list of all fastas in directory
    dbFileL = []
    for dbFile in glob.glob(fastaFilePath):
        dbStem = os.path.split(dbFile)[-1] # strain + _prot.fa or _dna.fa
        if "_prot.fa" in dbStem:
            dbStem = dbStem.split("_prot.fa")[0]
            if dbStem in strainNamesT:
                # this file name is in the tree, include it
                dbFileL.append(dbFile)
    return dbFileL
    
def formatDb(dbFileL,blastExecutDirPath):
    '''Format fasta files in fastaDir for blast by calling the formatdb
executable.'''
    for dbFileName in dbFileL:
        makeblastdbExecutable = os.path.join(blastExecutDirPath,'makeblastdb')
        subprocess.call([makeblastdbExecutable, '-dbtype' ,'prot', '-in', dbFileName],stdout=subprocess.PIPE)
    return

def makeBlastClineList(dbFileL_1, dbFileL_2,paramD):
    '''Create a list of lists, where the sublists have the command line
needed to run blastp on a pair of databases.'''

    fastaFilePath = paramD['fastaFilePath']
    blastFilePath = paramD['blastFilePath']
    blastExecutDirPath = paramD['blastExecutDirPath']

    # catch blast clines from old params files which don't have
    # trailing whitespace (can get rid of this eventually).
    if paramD['blastCLine'][-1] != ' ':
        blastCLine = paramD['blastCLine'] + ' ' + str(paramD['evalueThresh'])
    else:
        blastCLine = paramD['blastCLine'] + str(paramD['evalueThresh'])
    
    # get a tuple of the blastp command line args for use below
    blastCLineL = blastCLine.split()
    blastCLineL[0] = os.path.join(blastExecutDirPath,blastCLineL[0])
    blastCLineT = tuple(blastCLineL)

    # get blast file dir and extension
    splitT = os.path.split(blastFilePath)
    blastDir,rest = splitT
    blastExtension = rest.split("*")[-1]

    # get list of blast and fasta files currently
    fastaFiles = glob.glob(fastaFilePath)
    blastFiles = glob.glob(blastFilePath)

    # check if fasta files modified more recently than blast files
    shouldBlast = determineShouldBlast(dbFileL_1[0],dbFileL_2[0],blastFilePath)

    clineL=[]
    for query in dbFileL_1:
        for db in dbFileL_2:
            
            qstem = os.path.split(query)[-1]
            qstem = qstem.split("_prot.fa")[0]
            
            dbstem = os.path.split(db)[-1]
            dbstem = dbstem.split("_prot.fa")[0]
                
            outFN = os.path.join( blastDir, qstem + '_-VS-_' + dbstem + blastExtension )
                
            # only add to list if this blast file doesn't already exist
            if not os.path.isfile(outFN) or shouldBlast:
                L = list(blastCLineT) + ['-query',query,'-db',db,'-out',outFN]
                clineL.append(L)
                    
    return clineL

def determineShouldBlast(dbFile1,dbFile2,blastFilePath):
    '''Compare a blast and a fasta file. If the fasta file was edited
    more recently, return True, otherwise False.
    '''

    dbStem1 = os.path.split(dbFile1)[-1]
    dbStem1 = dbStem1.split("_prot.fa")[0] 

    dbStem2 = os.path.split(dbFile2)[-1]
    dbStem2 = dbStem2.split("_prot.fa")[0] 

    # get dirs and extensions
    blastDir,rest = os.path.split(blastFilePath)
    blastExtension = rest.split("*")[-1]

    fastaFN = dbFile1
    blastFN = os.path.join(blastDir, dbStem1 + '_-VS-_' + dbStem2 + blastExtension )

    if not os.path.isfile(blastFN):
        # hasn't been run, return True
        return True
    else:
        # times
        fastaModTime = os.path.getmtime(fastaFN)
        blastModTime = os.path.getmtime(blastFN)

        if fastaModTime < blastModTime:
            return False
        else:
            return True

def subprocessWrapper(cline):
    '''A wraper within which we call subprocess. Becasue we're having
blast write to file, we dump std_out. return std_err.'''
    pipes=subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = pipes.communicate()
    return stderr
