# Zunyan Wang
import sys
orthoFN=sys.argv[1]
blastFileList=sys.argv[2:]

def loadOrthos(orthoFN):
    '''Reads the orthologs file. The first lines of the file contains
species names, in the form: number species name. Gets all the names in
a list. Then reads the remaining lines to get ortholog sets. Puts each
set in a tuple, collecting the tuples in a list. Returns these two
lists.

    '''
    f= open(orthoFN,'r')

    strainL=[]
    orthoL=[]
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.rstrip().split()
        if s[0].isdigit(): # is a number-species line
            strainL.append(L[1])
        else:
            # is an ortho list line
            orthoL.append(tuple(L))
            
            
    f.close()
    return orthoL

def allPairs(orthoL):
    '''Get all pairs of genes in orthoL (lower triangular).'''
    groupOfPairs=[] 
    for orthoSet in orthoL: #loops over all sets of orthologs
    	for i in range(len(orthoSet)): #loops over every ortholog in a set
	    for j in range(len(orthoSet)): #pairing the previous ortholog with all other orthologs in this set
	    	if i<j:
		   groupOfPairs.append((orthoSet[i],orthoSet[j]))
                   groupOfPairs.append((orthoSet[j],orthoSet[i]))
    return set(groupOfPairs) #turning the result from a list to a set, which can speed up searching process.

def takeOutCoreScore(blastFileList, orthoL):
    '''Looks through a list of blast files and print out scores of similarity for the pairs of orthologs that we want.'''
    pairs=allPairs(orthoL) #get all the pairs that we want
    for blastFile in blastFileList: #loops through all blast files
        f=open(blastFile,'r')
        while True:
            s=f.readline()
            if s=='': #stop when reached the end of a file
                break
            L=s.rstrip().split("\t") #split a line into elements of a list
            if len(L)==12: #only keeps the line with actual data, which all have a length of 12
                firstStrain=L[0]
                secondStrain=L[1]
                if (firstStrain,secondStrain) in pairs:
                    print L[2] #print the similarity score if the pair of strain is what we look for
        f.close()
    return
	

takeOutCoreScore(blastFileList,loadOrthos(orthoFN))










    
