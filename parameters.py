# module for loading parameters.


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
