# Zunyan Wang
import sys

File = sys.argv[1]
numSpecies = int(sys.argv[2])
toPrint = int(sys.argv[2])

linesPerBlock = 2*numSpecies+1 # +1 for the trailing blank line

f= open(File, "r")

# This function prints out a specific block given an index out of all
# blocks of strains
def printBlock(index,blocks):
    rep=index*linesPerBlock
    for i in range(rep):
        blocks.readline()
    for j in range(linesPerBlock):
        line=blocks.readline()
        print line,
    return 

printBlock(toPrint,f)
f.close()
