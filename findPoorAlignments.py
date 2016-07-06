# Zunyan Wang
import sys, parseFuncs

File = sys.argv[1]

vF=int(sys.argv[2]) #Viewing Frame - number of bps looked at
mT=float(sys.argv[3]) #Mismatching threshold - % of the viewing frame must by the different
#The above two variables are parameters that user input. To get a comprehensive result, try a combination of viewing frames from 20 to 100, mismatching threshold from .4 to .5.

totalBlocks=801 #number of strains for every block
orthologBlockf = open(File, "r")
blockCount = 0 #keep track of which block currently examining
poorL=[] #list to store result

# Find poorly aligned regions in an ortholog block

def thresholdChecker(orthologBlock):
    """Checks if the ortholog block has regions with a matching percentage higher than the mT.  Will print the indices of the blocks that do are poor."""
    global blockCount
    for b in range(totalBlocks):
        seqL = parseFuncs.blockReader(orthologBlock) #blockReader is a helper function that converts sequence from .afa format into list format that we can program on.
        if seqL==[]: #check if sequence is empty
            print "bad"
            return False
        guideL = [] #store the result of comparing at each site into a list of characters, where 0 represent complete match, 1 represent mismatch, and 2 represent a dash(gap) is present at that site on any one strain.
        for k in range(len(seqL[0][1])): #in this loop, we are comparing the first strain to every other strain, one by one.
            match = True #default is complete match
            if seqL[0][1][k]=="-": #check whether first strain has a dash
                match="gap"
            else:
                for j in range(1, len(seqL)): #loop over every other strain
                    if seqL[0][1][k] != seqL[j][1][k]:
                        if seqL[j][1][k] == "-": #check whether comparing strain has a dash
                            match = "gap"
                            break
                        else: match = False #mismatch
            if match==True:
                guideL += [0] #storing result into guideL.
            elif match==False:
                guideL += [1]
            else:
                guideL += [2]
        #next section makes a moving window, recording the indices of blocks that contain windows with too much mismatches.
        good = True
        for i in [0,len(guideL),5] : #window moves at a rate of 5 sites per move
            seq=guideL[i:i+vF]
            countGap= count(seq, 2) #count is a helper function that count the number of a character, in this case 1 or 2.
            countMis= count(seq, 1)
            if countMis > mT*vF: #check if the number of mismatch is more than mT percent of the moving frame
                good=False
        if good == False:
            poorL.append(blockCount) #record indices of poor blocks
        blockCount+=1
    print poorL
    return poorL


def count(List, sym):
    """counts the number of symbols, sym,  present in the string Str"""
    counter = 0
    for i in range(len(List)):
        if List[i]==sym:
            counter += 1.0
    return counter

        
thresholdChecker(orthologBlockf)









    
        
