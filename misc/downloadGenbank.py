from Bio.Entrez.Parser import StringElement
from pathlib import Path
from Bio import Entrez
import sys, glob, os, ftplib, re
import warnings

xenoGI_path = os.path.join(sys.path[0],'..')
if xenoGI_path not in sys.path:
    sys.path.insert(0,os.path.join(sys.path[0],'..'))

from xenoGI import parameters

#---------------------------
#
# various variables
#
#---------------------------

# set global variables of paramFN and paramD
paramFN = 'params.py'
paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)

# set the global variable of the genbank directory
GBFFDIR:str = paramD['genbankFilePath'].split('/')[0]
# include the slash at the end of the genbank directory, so it's a folder
if GBFFDIR[-1]!='/':
    GBFFDIR += '/'

# set the global variable for the file map name
MAPFILENAME:str = paramD['fileNameMapFN']

# set the warnings format
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
        return '%s:%s: %s:%s\n' % (filename, lineno, category.__name__, message)

warnings.formatwarning = warning_on_one_line

#---------------------------
#
# from ncbiPythonTools
#
#---------------------------

def checkEntrezEmail(email:str) -> None:
    """ checkEntrezEmail:
            Accepts an email address as input. Checks if Entrez.email has been
            set. If not, then it attempts to set Entrez.email to the provided
            email address. If not possible, then it raises an exception. Does 
            not return.
    """
    # make sure Entrez.email can be set
    if Entrez.email is None:
        if email is not None:
            Entrez.email = email
        else:
            raise BaseException("Entrez.email has not been set.")

# Downloading GBFF

def downloadGbff(ftpPath:str, gbffDir:str) -> str:
    """ downloadGbff:
            Accepts two strings as inputs: one indicating the ftp path of the
            assembly and one indicating the directory where the downloaded file
            should be saved. Returns the filename of the newly downloaded file
            as a string.
    """
    # constants
    NCBI_FTP  = "ftp.ncbi.nlm.nih.gov"
    LEAD_STR  = "ftp://" + NCBI_FTP
    GREP_FIND = r'^.+/([^/]+)$'
    GREP_REPL = r'\1'

    # make sure gbff directory ends in a slash
    if gbffDir[-1] != '/':
        gbffDir += '/'
    
    # make the folder
    Path(gbffDir).mkdir(parents=True, exist_ok=True)

    # remove the leading string from the ftp path
    ftpPath = ftpPath[len(LEAD_STR):]

    # extract the filename from the ftp path
    gzFileName = re.sub(GREP_FIND, GREP_REPL, ftpPath)

    # add the folder info to the path
    gzFileName = gbffDir + gzFileName

    # download the file
    downloadFileFromFTP(NCBI_FTP, ftpPath, gzFileName)

    # get the filename for the unpacked gbff
    gbffFileName = removeFileExtension(gzFileName)

    # unpack the gz
    decompressGZ(gzFileName)

    # return the gbff filename (without directory information)
    return re.sub(GREP_FIND, GREP_REPL, gbffFileName)


def downloadFileFromFTP(ftpHost:str, ftpFilePath:str, filename:str) -> None:
    """ downloadFileFromFTP:
            Accepts three strings: the ftp server, the full path to the file to
            be downloaded, and the name of the file to write locally. Downloads
            the file and writes it. Does not return.
    """
    # initialize the ftp object, then connect and login to the server
    ftp = ftplib.FTP()
    ftp.connect(ftpHost)
    ftp.login()

    # construct the ftp retrieve command
    retrCmd = "RETR " + ftpFilePath

    # download the file and write it locally
    with open(filename, 'wb') as fileHandle:
        ftp.retrbinary(retrCmd, fileHandle.write)


def removeFileExtension(inFile:str) -> str:
    # remove the last period to the end of the string
    return re.sub(r'\.[^\.]+$', '', inFile)


def decompressGZ(gzFileName:str) -> None:
    """ decompressGZ:
            Accepts a compressed gbff file as input. Extracts the gbff file and
            deletes the compressed file. Does not return.
    """
    # command construction
    cmd = 'gzip -d ' + gzFileName
    
    os.system(cmd)

#---------------------------
#
# helper functions
#
#---------------------------

def ftpToFileName(ftpPath:str):
    '''takes in an ftp path, in the format given in the Entrez ncbi assembly 
       summary, and returns the ftp path joined to the file name of the genbank
       of the genome. This filename can be taken in by downloadGbff'''
    # only take in strings
    if type(ftpPath) != str and type(ftpPath) != StringElement:
        raise TypeError("ftpToFileName can only take in a string")
    # get just the folder where the file is stored
    folderName = ftpPath.split('/')[-1]
    # add the proper ending to get the full genbank file
    fileName = folderName + '_genomic.gbff.gz'
    # return the full path to the file
    return ftpPath + '/' + fileName

def nonAssemblyIdToFTP(idNum, database:str):
    '''takes in id(s) (either a single string for one id number or a list of id
       numbers) and database, and returns the ftpPaths and species names for each
       id. The ids must be given in as strings, not ints'''
    # criteria for valid assembly links
    filter = '"representative genome"[filter] OR "reference genome"[filter]'
    # get an elink from the idNum to assembly db
    handle = Entrez.elink(dbFrom=database, id = idNum, db='assembly', term = filter)
    linkResults = Entrez.read(handle, validate = False)
    handle.close()
    # initialize list for storing assembly numbers
    assemblyIds = []
    # initialize list for retries
    retryIds = []
    # iterate through results
    for i in range(len(linkResults)):
        # pull an assembly number for each result
        try:
            assemblyIds.append(linkResults[i]['LinkSetDb'][0]['Link'][0]['Id'])
        except:
            retryIds.append(linkResults[i]['IdList'])
    # retry ids without filter
    if retryIds != []:
        handle = Entrez.elink(dbFrom=database, id = retryIds, db='assembly')
        linkResults = Entrez.read(handle, validate = False)
        handle.close()
        for i in range(len(linkResults)):
            # pull an assembly number for each result
            try:
                assemblyIds.append(linkResults[i]['LinkSetDb'][0]['Link'][0]['Id'])
            except:
                raise Exception("No assembly entries were found for the id " + str(linkResults[i]['IdList']))
    # give assembly numbers to assemblyIdToFTP
    return assemblyIdToFTP(assemblyIds)

def genomeIdToFTP(idNum, database:str):
    '''takes in genome id(s) (either a single string for one number or a list 
       of id numbers), and returns the ftpPaths as a list, and their corresponding
       species names as a list. The ids must be given in string form, not as ints'''
    # initialize a list for results
    dbEsummary = []
    # sort input types
    if type(idNum) == str:
        handle = Entrez.esummary(id = idNum, db=database)
        dbEsummary = Entrez.read(handle, validate = False)
        handle.close()
    elif type(idNum) == list:
        # possible errors with multiple inputs being given in as python list
        # change python list into a single comma delimited string
        stringAssemblyID = ",".join(idNum)
        # give string to query
        handle = Entrez.esummary(id = stringAssemblyID, db=database)
        dbEsummary = Entrez.read(handle, validate = False)
        handle.close()
    # initialize list for accession numbers
    assemblyIds = []
    # go through each entry, and get the assembly id
    for dict in dbEsummary:
        # add the assembly id to the list
        assemblyIds.append(dict['AssemblyID'])
    # get ftp's and species names from assembly id's
    return assemblyIdToFTP(', '.join(assemblyIds))

def assemblyIdToFTP(assemblyID):
    '''takes in assembly id(s) (either a single string for one number or a list 
       of id numbers), and returns the ftpPaths as a list, and their corresponding
       species names as a list. The ids must be given in string form, not as ints'''
    # get esummaries for each assembly id
    try:
        # directly paste assembly ids into query
        handle = Entrez.esummary(id = assemblyID, db='assembly')
        assemblyEsummary = Entrez.read(handle, validate = False)
        handle.close()
    except:
        # possible errors with multiple inputs being given in as python list
        # change python list into a single comma delimited string
        stringAssemblyID = ", ".join(assemblyID)
        # give string to query
        handle = Entrez.esummary(id = stringAssemblyID, db='assembly')
        assemblyEsummary = Entrez.read(handle, validate = False)
        handle.close()
    # initialize lists for ftps and species names
    ftpList = []
    speciesList = []
    # get the number of results from the esummary
        # note: should be the same as number of entries, may not be
    numResults = len(assemblyEsummary['DocumentSummarySet']['DocumentSummary'])
    # if numResults is not the number of id's given in, raise an exception
    if len(assemblyID) != numResults:
        checkMissing(assemblyID, assemblyEsummary)
    # iterate through the results
    for i in range(numResults):
        # get the current entry
        currentEntry = assemblyEsummary['DocumentSummarySet']['DocumentSummary'][i]
        # check that the entry is valid
        validateSummary(currentEntry)
        # get the ftp and species of the result
        currentFTP = currentEntry['FtpPath_GenBank']
        currentSpecies = currentEntry['SpeciesName']
        if len(currentEntry['Biosource']['InfraspeciesList']) != 0:
            if 'Sub_value' in currentEntry['Biosource']['InfraspeciesList'][0].keys():
                currentSpecies += '_' + currentEntry['Biosource']['InfraspeciesList'][0]['Sub_value']
        # replace spaces with underscores and dashes with underscores for current species
        currentSpecies = currentSpecies.replace(' ', '_')
        currentSpecies = currentSpecies.replace('-', '_')
        # do not add the ftp to the list of ftps if it is already in the list
        if currentFTP in ftpList:
            warnings.warn('the path for ' + currentSpecies + ' is a duplicate: skipping')
        else:
            # add the ftp
            ftpList.append(currentFTP)
            # check that the species is not already in speciesList
            if currentSpecies in speciesList:
                # if the species name is a duplicate, get a new file name
                numberedName = newFileName(currentSpecies, speciesList)
                # replace old name with better, numbered name
                currentSpecies = numberedName
            # add the species name
            speciesList.append(currentSpecies)
    
    # inbuilt test
    if len(ftpList) != len(speciesList):
        raise Exception('length of ftpList and speciesList are unequal')
    
    # return all ftps and species
    return (ftpList, speciesList)

def checkMissing(inputIds, summary):
    '''checks a string of ids against a related assembly summary, and issues a 
       warning of which ids do not have an associated summary entry. This function
       only works with assembly summaries'''
    # if input is string, convert to list
    if type(inputIds) == str:
        inputIds = inputIds.split(', ')
    # make a copy of inputIds, so it will not be edited
    missingIds = inputIds
    # get all the id numbers in the summary
    for i in range(len(summary['DocumentSummarySet'])):
        # get the id number for this entry
        idNum = summary['DocumentSummarySet']['DocumentSummary'][i].attributes['uid']
        # remove the idnum from the list of input ids
        missingIds.remove(idNum)
    if len(missingIds) > 0:
        # what's left in the input ids is what was not found. Convert to string
        stringIds = ', '.join(missingIds)
        warnings.warn('Assembly entries were not found for the numbers '+stringIds)

def newFileName(oldName:str, listOfNames:list) -> str:
    '''Takes in a string/name, and a list of names. If the given name already
       exists in the list, an underscore and a number are addeed to the end of the 
       list. This name is also checked for list uniqueness. The number increases
       each time through the loop, until we have a new name that is unique'''
    # start with the old name
    proposedName = oldName
    # pick a number to add to the end
    num = 1
    # create a loop for finding a new name, break the loop when the name is unique
    while proposedName in listOfNames:
        # create a new name with form oldName_num
        proposedName = oldName + '_' + str(num)
        # if the name is not unique, increment num and get a new name
        num += 1
    return proposedName

def validateSummary(summary):
    '''Takes in the entry for a single summary, which is given by 
       esummary['DocumentSummarySet']['DocumentSummary'][int] . This function will
       raise an exception if the entry does not have annotations or is not a 
       full-genome-representation, and will issue a warning if it does not have
       "Complete Genome" as its status'''
    # get the id number for communication
    idNum:str = str(summary.attributes['uid'])
    # check for annotation
    if 'has_annotation' not in summary['PropertyList']:
        # raise an exception
        raise Exception("The genome for id number " + idNum + " has no annotations")
    # check for full genome
    if 'full-genome-representation' not in summary['PropertyList']:
        # raise an exception
        raise Exception("The genome for id number " +idNum+ " is not a full genome.")
    # check for assembly status
    if summary['AssemblyStatus'] != 'Complete Genome':
        # issue a warning
        warnings.warn("The genome for id number " +idNum+ " is not a 'Complete Genome'."\
            + " Downloading it anyways")

def convertInputsToUIDs(idNums:list, database:str, email=None):
    '''Given a list of id numbers (uid's OR accession numbers, or a mix of) returns
       a list of all the idnumbers as UID's (converts accession numbers to UIDS)'''
    # initialize a string to use in Entrez search
    uids = ''
    # iterator through id's
    for id in idNums:
        # check if it is a digit uid
        if id.isdigit():
            # add as a uid
            uids += 'OR ' + id + '[uid] '
        else:
            # add as all fields
            uids += 'OR ' + id + '[All Fields] '
    # set Entrez email
    if email != None:
        Entrez.email = email
    # run search, changing retmax to be the number of inputs
    handle = Entrez.esearch(term = uids, db = database, retmax = len(idNums))
    result = Entrez.read(handle, validate = False)
    handle.close()
    # get just the ids
    idList = result['IdList']
    # check all were found
    if len(idList) != len(idNums):
        warnings.warn(""" The number of results found was not the number of id's given.
        Id's given: """ + str(idNums) + """
        Id's found: """ + str(idList))
    return idList

def getInfoFromAssemblyUID(idNums:list, email = None):
    '''Takes in a list of UIDs only, gets their associated assembly entries, and 
       returns the accession numbers and species names for each entry'''
    # initialize place for accession numbers
    accNums = []
    # initialize place for species/strain name
    speciesNames = []
    # check email
    if email != None:
        checkEntrezEmail(email)
    # pull a summary for each id
    handle = Entrez.esummary(id = ', '.join(idNums), db = 'assembly')
    Esummary = Entrez.read(handle, validate = False)
    handle.close()
    # get the accNum from each summary
    for documentSummary in Esummary['DocumentSummarySet']['DocumentSummary']:
        # add the accession number to the list
        accNums.append(documentSummary['AssemblyAccession'])
        currentSpecies = documentSummary['SpeciesName']
        if len(documentSummary['Biosource']['InfraspeciesList']) != 0:
            if 'Sub_value' in documentSummary['Biosource']['InfraspeciesList'][0].keys():
                currentSpecies += '_' + documentSummary['Biosource']['InfraspeciesList'][0]['Sub_value']
        # replace spaces with underscores and dashes with underscores for current species
        currentSpecies = currentSpecies.replace(' ', '_')
        currentSpecies = currentSpecies.replace('-', '_')
        speciesNames.append(currentSpecies)
    return accNums, speciesNames

def checkOutgroup(speciesList):
    '''In params.py there should be an outgroup listed; checkOutgroup searches 
       the given speciesList for that outgroup. If the outgroup is not found, 
       asks for userInput. If no outgroup is set, again asks for input (empty,
       does nothing to ensure user makes active choice)'''
    # attempt to get the outgroup
    try:
        # pull the outgroup
        outgroup = paramD['outGroup']
        if outgroup == '':
            raise KeyError
        # if the outgroup is not found
        if outgroup not in speciesList:
            # outgroup not found message
            print('\nThe outgroup ' +outgroup+ ' as provided in params.py does not appear in the results.' +
                  'This may be because the name does not match exactly, or because no analog exists.')
            userChoice = input('Would you like the outgrouop to be searched for and have it downloaded now? [y/N] ')
            # search for outgroup
            if userChoice.lower() == 'yes' or userChoice.lower() == 'y':
                # find the outgroup
                id = searchAssembly(keyword=outgroup, retmax = 1)
                # get the name from the outgroup
                ftp, specName = assemblyIdToFTP(id)
                # get the downloadable file name from the ftp path
                fixedFtp = ftpToFileName(ftp[0])
                # get the correct out file name
                outFileName = GBFFDIR + removeFileExtension(fixedFtp.split('/')[-1])
                # check that the file doesn't exist before downloading
                if glob.glob(outFileName) == []:
                    # download file
                    downloadGbff(fixedFtp, GBFFDIR)
                    # append to human map file
                    makeHumanMap([outFileName], [outgroup], MAPFILENAME, 'a')
                # file already exists
                else:
                    print('\nthis file already exists: please change either the name' +
                          ' in params.py or the human map file. As a reminder, this' +
                          ' name should contain no dashes, spaces, commas or special' +
                          ' characters')
            # user does not wish to download outgroup
            else:
                print('\nPlease either manually download and add the outgroup from' +
                      ' params.py to the human map file, or choose one of the already' +
                      ' downloaded files to be the outgroup - setting the outgroup' +
                      ' variable in params.py to be the human-readable name in the' +
                      ' human map file')
    except KeyError:
        input('\nNo outgroup has been specified in the given params.py. If no species' +
              ' tree is provided, xenoGI will attempt to make one and will require' +
              ' an outgroup to do so. Please either choose one of the downloaded files' +
              ' to be the outgroup - setting the outgroup variable in params.py to be' +
              ' the human-readable name in the human map file, or provide a tree. ')

#---------------------------
#
# downloads
#
#---------------------------

def downloadMultipleGBFFs(listOfIdNums:list, database:str, email:str = None, listOfHumanNames:list = None):
    '''Takes in a list of id's and their database and downloads all
       files. A map file (file names to human readable names) is made automatically
       using the species names as defined in the assembly entrez esummary. Optionally,
       a list of names to be used instead of the species names can be specified as
       listOfHumanNames. Likewise, an email is also optional, and if provided will
       be used to interact with Entrez. If no email is specified, Entrez will still
       run (albeit begrudgingly).'''
    # GET FTPS

    # set up entrez email
    checkEntrezEmail(email)
    # if we've got non-assembly numbers, run elink on them
    if database != 'assembly':
        ftps, specNames = nonAssemblyIdToFTP(listOfIdNums, database)
    # if we've got assembly numbers, go straight to summary
    else:
        ftps, specNames = assemblyIdToFTP(listOfIdNums)
    
    # DOWNLOAD FTPS

    # initialize file storage
    fileNames = []
    # iterate through the ftp path list
    for i in range(len(ftps)):
        # get the downloadable file name from the ftp path
        fixedFtp = ftpToFileName(ftps[i])
        # get the correct out file name
        outFileName = GBFFDIR + removeFileExtension(fixedFtp.split('/')[-1])
        # add the outfile name to the list of file names
        fileNames.append(outFileName)
        # check that the file doesn't exist before downloading
        if glob.glob(outFileName) == []:
            downloadGbff(fixedFtp, GBFFDIR)
        
        # if it already exists, we indicate it

        else:
            warnings.warn('file '+outFileName+' already exists. Skipping')

    # HUMAN MAP

    # if human names specified, we fill out any files w/out names
    if listOfHumanNames != None:
        '''iterate through, for each '' entry, replace with specnames'''
        for i in range(len(listOfHumanNames)):
            if listOfHumanNames[i] == '':
                listOfHumanNames[i] = specNames[i]
        specNames = listOfHumanNames
    # check if the mapfilename exists. If it doesn't, make it
    if glob.glob(MAPFILENAME) == []:
        makeHumanMap(fileNames, specNames, MAPFILENAME)
    # if the file does exist, we let the user know
    else:
        # ask the user if they want to write over it 
            # NOTE: we do this because makeHumanMap is write-only, not append
        userChoice = input('The filename-humanname mapping file already exists, would you like to write over it? [Y/n]  ') or 'yes'
        # if the user indicates yes, Yes, y, or Y, we overwrite the file
        if userChoice.lower() == 'yes' or userChoice.lower() == 'y':
            makeHumanMap(fileNames, specNames, MAPFILENAME)
        else:
            print('the file was not rewritten')
    
    # warn about outgroups!
    checkOutgroup(specNames)

#---------------------------
#
# wrapper functions
#   user input -> downloads
#
#---------------------------

def searchToDownload(handle, database, email):
    '''taking in the handle of an already-completed search, downloads the files
       using downloadMultipleGBFFS'''
    # give downloadMultipleGBFFs the IdList
    downloadMultipleGBFFs(handle['IdList'], database)

def queryToDownload(keyword:str, retmax:int=10, email = None):
    '''This function takes in a single string (the search) and searches ncbi's 
       assembly database for matches. This search is also filtered so all results
       have annotations and are a full genome representation. Additionally, the
       results are sorted by significance. The results are then downloaded as genbanks,
       and a human map file is made as well. retmax defaults to 10'''
    # define variable for output file name
    outputFN = 'searchResults.txt'
    # search assembly database
    ids = searchAssembly(keyword, retmax, email)
    # turn the ids into accessions
    accNums, specNames = getInfoFromAssemblyUID(ids)
    # write the accNums and specNames to a file
    makeHumanMap(accNums, specNames, outputFN)
    # get user input to check file
    userChoice = input('Search results were written to the file '+outputFN+' which is editable.' +
                       ' If there are any changes you would like to make, such as changing names, adding or deleting results,' +
                       ' or adding an outgroup, please do so now. Any changes made before continuing will be used when downloading.' +
                       ' Would you like to download the results as shown in the (edited or unedited) file? [y/N] ') or 'no'
    if userChoice.lower() == 'yes' or userChoice.lower() == 'y':
        fileToDownload(outputFN, 'assembly', email)
    else:
        print('the search results were not downloaded')

def searchAssembly(keyword:str, retmax:int = 10, email=None):
    '''This function takes in a single string (the search) and searches ncbi's 
       assembly database for matches. This search is also filtered so all results
       have annotations and are a full genome representation. Additionally, the
       results are sorted by significance. If no retmax is given, default is 10.
       The id numbers are returned'''
    # check email
    checkEntrezEmail(email)
    # input desired filter for annotations and full genome representation
    properties = ' AND "has annotation"[Properties] AND "full genome representation"[filter]'
    # create full search term
    searchTerm = keyword + properties
    # run the search
    handle = Entrez.esearch(db = 'assembly', term = searchTerm, retmax = retmax, sort = 'Significance')
    result = Entrez.read(handle, validate = False)
    handle.close()
    # get the ids from the search
    ids = result['IdList']
    return ids

def fileToDownload(mapFileName, database, email = None):
    '''takes in a map file with each line having expected format 
       idnum \t readablename. From there, downloads the corresponding genbank
       with each of the idnums as they correspond with the database, and make
       a humanMap file using the names provided. idnums that do not have names
       will have names provided for them. The email is for Entrez purposes.'''
    # initialize list of idnums and species names
    idNums = []
    speciesNames = []
    # read the file
    with open(mapFileName, 'r') as file:
        while True:
            # pop the first line
            oneLine = file.readline()
            # if empty file, break the loop
            if oneLine == '':
                break
            # change fileLine from string to list of strings
            args = oneLine.strip().split('\t')
            # rename for easier working, and strip it to get rid of any extra spaces
            currentNum = args[0].strip()
            # if no species is given, we add the number and use '' for the name
            if len(args) == 1:
                idNums.append(currentNum)
                speciesNames.append('')
            
            # if a species name is given, it must be handled.
            else:
                # rename for easier working, and strip it to get rid of outside spaces
                currentName = args[1].strip()
                # change inner spaces to underscores
                currentName = currentName.replace(' ', '_')
                # case for unique current name, or empty current name having been stripped
                if currentName not in speciesNames or currentName == '':
                    # add the id number and species number
                    idNums.append(currentNum)
                    speciesNames.append(currentName)
                # if the current name is not unique, we skip or rename
                else:
                    # ask the user if they want to skip or rename
                    userChoice = input('The human-readable name ' + currentName + ' already exists. ' +
                                       'Would you like to skip this entry, or have it renamed? [skip/rename] ')
                    validInput = False
                    # we can only move on if the response was valid
                    while validInput == False:
                        # user chooses skip
                        if userChoice.lower() == 'skip':
                            # break the loop! We are skipping
                            validInput = True
                        # user chooses a rename
                        elif userChoice.lower() == 'rename':
                            numberedName = newFileName(currentName, speciesNames)
                            # append the new species name
                            speciesNames.append(numberedName)
                            # append the current number
                            idNums.append(currentNum)
                            # break the loop
                            validInput = True
                        else:
                            # neither skip or rename, ask for the input again
                            userChoice = ('Response was not not understood: please reply either "skip" or "rename"')
    # check the entrez email
    checkEntrezEmail(email)
    # here, ID numbers could be accession or uid's: we run a search to get just uid's
    uids = convertInputsToUIDs(idNums, database, email)
    # use idNums and speciesNames to perform downloads and mapfile making
    downloadMultipleGBFFs(uids, database, listOfHumanNames=speciesNames, email=email)

def _summaryToSpeciesName(handle):
    '''given an already completed assembly esummary, return the first species name'''
    # pull the species name from the handle
    speciesName = handle['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
    return speciesName.replace(' ', '_')

def makeHumanMap(fileNames:list, speciesNames:list, outputName, type='w'):
    '''use the output from downloadMultipleGBFFs to make the map file relating
       file names and their species names'''
    # read in the file as 'w', creating it if it does not exist
    with open(outputName, type) as f:
        # iterate through the file names
        for i in range(len(fileNames)):
            # add 'fileName tab speciesName newline' to file
            f.write(fileNames[i].split('/')[-1] + "\t" + speciesNames[i] + "\n")
    # with block closes the file automatically

#------------------------
# 
# Main function
# (for cmd line arguments)
# 
#------------------------

def main():
    '''a main function to take in command line arguments, and run the download
       without entering the interactive mode. The first argument taken in is the
       file name of a file containing id numbers in the first column, and optional
       names in the second. The second argument is the database associated with 
       the id numbers, and should be in a form Entrez takes in (link below). The
       third, optional argument is an email for using Entrez. Entrez will run
       (begrudgingly) if no email is given, but the option is there if so desired.'''
    # ensure we have either two or zero arguments
    try:
        # check the number of arguments, where the first argument is always downloadGenbank
        assert(len(sys.argv) == 1 or len(sys.argv)==4)
    # if incorrect number of arguments given, quit
    except:
        warnings.warn(
            """

    Either three or zero parameters are required. If three, must be
        1. A text file of id numbers, and optionally their human-readable names
        2. The database the id numbers correspond to which must be one of the E-utility Database Names from
        https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
        3. An email to make Entrez calls with
    For example:
        downloadGenbank.py nucleotideIDfile.txt nuccore username@email.com
            """
            )#,file=sys.stderr)
        sys.exit(1)
    
    # three arguments -> file name, database, email
    if len(sys.argv) == 4:
        # read in the arguments
        fileName = sys.argv[1]
        database = sys.argv[2]
        email = sys.argv[3]
        # download, with email!
        fileToDownload(fileName, database, email)
    # no arguments -> set up interactive console
    else:
        # import required packages
        import readline, code, rlcompleter
        # code stolen from interactiveAnalysis xenoGI function
        vars = globals()
        vars.update(locals())
        readline.set_completer(rlcompleter.Completer(vars).complete)
        readline.parse_and_bind("tab: complete")
        code.InteractiveConsole(vars).interact()

# RUN MAIN !
if __name__ == "__main__" :
    main()
