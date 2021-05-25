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

    blastFilePath = paramD['blastFilePath']
    blastExecutDirPath = paramD['blastExecutDirPath']
    blastFileJoinStr = paramD['blastFileJoinStr']
    
    # catch blast clines from old params files which don't have
    # trailing whitespace (can get rid of this eventually).
    if paramD['blastCLine'][-1] != ' ':
        blastCLine = paramD['blastCLine'] + ' ' + str(paramD['evalueThresh'])
    else:
        blastCLine = paramD['blastCLine'] + str(paramD['evalueThresh'])

    # get a tuple of the blastp command line args for use below
    blastCLineL = processCline(blastCLine)
    blastCLineL[0] = os.path.join(blastExecutDirPath,blastCLineL[0])
    blastCLineT = tuple(blastCLineL)

    # get blast file dir and extension
    splitT = os.path.split(blastFilePath)
    blastDir,rest = splitT
    blastExtension = rest.split("*")[-1]

    # check if fasta files modified more recently than blast files
    shouldBlast = determineShouldBlast(dbFileL_1[0],dbFileL_2[0],blastFilePath,blastFileJoinStr)

    clineL=[]
    for query in dbFileL_1:
        for db in dbFileL_2:

            # split of paths and extensions
            qstem = os.path.split(query)[-1]
            if "_prot.fa" in qstem:
                qstem = qstem.split("_prot.fa")[0]
            else:
                # just remove extension
                qstem = os.path.splitext(qstem)[0]
                
            dbstem = os.path.split(db)[-1]
            if "_prot.fa" in dbstem:
                dbstem = dbstem.split("_prot.fa")[0]
            else:
                dbstem = os.path.splitext(dbstem)[0]
                
            outFN = os.path.join( blastDir, qstem + blastFileJoinStr + dbstem + blastExtension )
                
            # only add to list if this blast file doesn't already exist
            if not os.path.isfile(outFN) or shouldBlast:
                L = list(blastCLineT) + ['-query',query,'-db',db,'-out',outFN]
                clineL.append(L)
                    
    return clineL

def processCline(cline):
    '''Process the command line parameter so that the output format
section can be passed in to subprocess as a unit.
    '''
    # split out the part between quotes
    L = cline.split('"') # " make emacs happy.
    return L[0].split() + [L[1]] + L[2].split()
    
def determineShouldBlast(dbFile1,dbFile2,blastFilePath,blastFileJoinStr):
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
    blastFN = os.path.join(blastDir, dbStem1 + blastFileJoinStr + dbStem2 + blastExtension )

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

def createBlastD(strainNamesT,blastFileJoinStr,blastDir,blastExt,evalueThresh,alignCoverThresh,percIdentThresh):
    '''Create a dictionary of blast output. Keyed by strain pair, with
values being the list returned by parseBlastFile.
    '''
    pathL,strainPairL = constructBlastFileList(strainNamesT,blastFileJoinStr,blastDir,blastExt)

    # load all blast files
    blastD = {}
    for i in range(len(pathL)):
        blastPath = pathL[i]
        strainPair = strainPairL[i]
         # set coverage and percid thresh to 0, to take all.
        blastD[strainPair] = parseBlastFile(blastPath,evalueThresh,alignCoverThresh,percIdentThresh)

    return blastD,strainPairL

def constructBlastFileList(strainNamesT,blastFileJoinStr,blastDir,blastExt):
    '''Constructs the path to all blast files that should be present given
the strains in strainNamesT. Returns two lists of the same
length. pathL, with the paths. strainPairL with tuples of the pairs of
strains.

    '''
    pathL = []
    strainPairL = []
    for strain1 in strainNamesT:
        for strain2 in strainNamesT:
            path = os.path.join(blastDir,strain1+blastFileJoinStr+strain2+blastExt)
            pathL.append(path)
            strainPairL.append((strain1,strain2))
    return pathL,strainPairL

def parseBlastFile(blastFN,evalueThresh,alignCoverThresh,percIdentThresh):
    '''Parse a single blast output file, returning all hits as a list of
tuples. alignCoverThresh is a threshold for the length of the
alignment relative to query and subject length. The files have the
following fields (our custom output) qseqid sseqid evalue qlen qstart
qend slen sstart send pident score.
    '''
    L=[]
    with open(blastFN,'r') as f:
        s = f.readline()
        while s:
            blastLine = s.split()

            # get numeric gene from front, if not present get whole name
            queryStr = blastLine[0].split('_')[0]
            if queryStr.isdigit():
                queryGene = int(queryStr)
            else:
                queryGene = blastLine[0]

            subjectStr = blastLine[1].split('_')[0]
            if subjectStr.isdigit():
                subjectGene = int(subjectStr)
            else:
                subjectGene = blastLine[1]

            evalue = float(blastLine[2])
            if evalue < evalueThresh:
                qlen = int(blastLine[3])
                qstart = int(blastLine[4])
                qend = int(blastLine[5])
                slen = int(blastLine[6])
                sstart = int(blastLine[7])
                send = int(blastLine[8])
                pident = float(blastLine[9])
                score = int(blastLine[10])

                alCov = ((qend-qstart) + (send-sstart)) / (qlen+slen)
                alLen = ((qend-qstart) + (send-sstart)) / 2 # av align len
                
                if alCov > alignCoverThresh and pident > percIdentThresh:
                    L.append((queryGene,subjectGene,evalue,alCov,pident,score,alLen))
            s = f.readline()
    return L

def getBestHitsDictionary(blastL):
    '''getBestHitsDictionary: Accepts a parsed blast table as input (list
            of tuples). Uses the 'score' to determine the best hit for
            each unique query. Returns a dictionary containing the
            best hits with the query names as the keys and
            (hit,pident,score) as value.
    '''
    # constants
    QUERY_INDEX = 0
    HIT_INDEX = 1
    PIDENT_INDEX = 4
    SCORE_INDEX = 5
    ALIGNLEN_INDEX = 6
    
    # initialize the output
    bestHitD = {}

    for entry in blastL:
        # put values into more compact variable names
        query = entry[QUERY_INDEX]
        hit = entry[HIT_INDEX]
        pident = entry[PIDENT_INDEX]
        score = entry[SCORE_INDEX]
        alLen = entry[ALIGNLEN_INDEX]
        
        # only one entry per query, only keep the one with the highest score
        if query not in bestHitD.keys() or score > bestHitD[query][2]:
            bestHitD[query] = (hit,pident,score,alLen)
    
    return bestHitD

