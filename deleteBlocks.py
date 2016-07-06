# Zunyan Wang
import sys

File = sys.argv[1]
numBlocks=int(sys.argv[2])
numSpecies = int(sys.argv[3])
toDelete = list(sys.argv[4:])
toDeleteList=map(int,toDelete)

linesPerBlock = 2*numSpecies+1 # +1 for the trailing blank line

f= open(File, "r")

def deleteBlocks(indices,numBlocks,blocks):
    for i in range(numBlocks):
        if i in indices:
            for j in range(linesPerBlock):
                line=blocks.readline()
        else:
            for k in range(linesPerBlock):
                line=blocks.readline()
                print line,
    return

deleteBlocks(toDeleteList, numBlocks, f)
f.close()
