# module for loading parameters.
import os

def loadParametersD(paramFN):
    '''Load a parameters file into a dictionary. Assumes the file consists
of python assignment statements.'''
    
    paramD={}
    f=open(paramFN,'r')

    while True:
        s=f.readline()
        if s == '':
            break
        if s[0]=='#' or s[0] == '\n':
            continue
        # if it's not the EOF, a comment, or a blank line, so it
        # must be an assignment statement.
        key,value = s.rstrip().split('=')
        key = key.strip()
        value = value.strip()
        paramD[key] = eval(value)

    f.close()

    return paramD

def loadFileNameMapD(fileNameMapFN,genbankFileList=None):
    '''Create a dictionary with mappings between genbank file names and
the human readable names we use in the tree. If fileNameMapFN contains
a string, we load that file and construct the mappings based on this.
Expects file with one species per line: genbank name + white space +
human name. If fileNameMapFN is None, we create the mappings between
the full file names (with path and extension) and the stem of the file
name, which in this case should correspond to what is in the input
tree.
    '''
    fileNameMapD = {}
    if fileNameMapFN == None:
        for fullPathFN in genbankFileList:
            fn = os.path.split(fullPathFN)[-1]
            stem = os.path.splitext(fn)[0]
            fileNameMapD[fn] = stem
    else:
        f = open(fileNameMapFN,'r')
        while True:
            s = f.readline()
            if s == '':
                break
            elif s[0].isspace():
                continue
            genbankStem,human = s.rstrip().split()
            fileNameMapD[genbankStem] = human
    return fileNameMapD
