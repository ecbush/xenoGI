import sys,os
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,genomes,trees,families,scores,islands
from analysis import *

def islandsOfInterest(minLength):
    longIslands = islandsInStrainLongEnough(minLength)
    longOnChrom = islandsOnChromosome(longIslands)
    longOnChromInRange,overlapList,totalBases,islandsList = islandsInRange(longOnChrom)
    overlap = sum(overlapList)
    print("Islands: ", islandsList)
    print("overlap: ", overlapList)
    print("total bases: ", totalBases)
    print('pecent overlap: ', overlap/totalBases)
    return longOnChromInRange

def islandsInStrainLongEnough(minLength):
    '''returns a list of the xenoGI islands that are longer than minLength'''
    
    #list of islands in strain, empty list of islands to return
    potentialIslands = islandByNodeL[node]
    returnIslands = []
    
    #for each island, check that it meets our criteria
    for island in potentialIslands:
        islandGenesInStrainL = getIslandGenesInStrain(island,strainNum,familyT)
        
        #check island's length > minLength
        neighbGenesL,firstIslandGene,lastIslandGene=getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0)
        startPos =min(int(geneInfoD[geneNames.numToName(firstIslandGene)][4]), int(geneInfoD[geneNames.numToName(firstIslandGene)][5]))
        endPos = max(int(geneInfoD[geneNames.numToName(lastIslandGene)][5]),int(geneInfoD[geneNames.numToName(lastIslandGene)][4]))
        #if the island is long enough, add it to our list of potential islands
        if (endPos-startPos)>minLength: returnIslands.append(island)
    return returnIslands
    
def islandsOnChromosome(potentialIslands):
    returnIslands = []
    #loop through each island,if an island is on the correct chromosome,
    #add it to our list of potential islands
    for island in potentialIslands:
        islandGenesInStrainL = getIslandGenesInStrain(island,strainNum,familyT)
        chromFound = geneInfoD[geneNames.numToName(islandGenesInStrainL[0])][3]
        if (chromFound == chrom): returnIslands.append(island)
    return returnIslands

def islandsInRange(potentialIslands):
    validationRanges, totalBases=readRanges(validationFile)
    returnIslands = []  #holds islands that have any overlap
    overlapList=[] #holds the # of bp that overlap w/ a validation island for each island we check
    islandsList=[0]*len(validationRanges) #holds the number of xenoGI islands that overlap w/each validation range
    
    for island in potentialIslands:
        #get the start and end position for the islands
        islandGenesInStrainL = getIslandGenesInStrain(island,strainNum,familyT)
        neighbGenesL,firstIslandGene,lastIslandGene=getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0)
        startPos = min(int(geneInfoD[geneNames.numToName(firstIslandGene)][4]), int(geneInfoD[geneNames.numToName(firstIslandGene)][5]))
        endPos = max(int(geneInfoD[geneNames.numToName(lastIslandGene)][5]),int(geneInfoD[geneNames.numToName(lastIslandGene)][4]))

        #check that island is in validation range
        inRange,overlap,indices=islandInRange(validationRanges, startPos, endPos)

        #if the island overlaps with any validation range, print an array showing which ones
        #this will be used to see which islands cover multiple ranges
        if sum(indices)>1: print(island,indices)
        #update overlap list
        overlapList.append(sum(overlap))
        #update the islandsList so it reflects how many xenoGI islands are in each range
        islandsList=list(map(lambda x,y:x+y, islandsList, indices))
        #add the island to the returnIslands list if it overlaps with any validation range
        if inRange:
            returnIslands.append(island)

    return returnIslands, overlapList, totalBases, islandsList

def islandInRange(validationRanges, startPos, endPos):
    '''returns true if the island overlaps with any 
    validation ranges and false otherwise'''
    overlap = [0]*len(validationRanges)
    indices = [0]*len(validationRanges)
    for index in range(0,len(validationRanges)):
        vRange=validationRanges[index]
        overlap[index]=min((vRange[1]-vRange[0]),max((endPos-vRange[0]),0))-min((vRange[1]-vRange[0]),max((startPos-vRange[0]),0))
        if overlap[index]>0: indices[index]=1
    if sum(overlap)>0: return True, overlap, indices
    return False, overlap, indices

def readRanges(fileName):
    '''Read a file with start end ranges separated by tabs, and return list of lists.'''
    totalBases = 0
    f=open(fileName,'r')
    outLL = []
    while True:
        s=f.readline()
        if s == "":
            break
        L=[]
        for elem in s.split('\t'):
            L.append(int(elem))
        totalBases += (L[1]-L[0])
        outLL.append(L)
    return outLL, totalBases

def nodesPerRange(fileName):
    '''read a tab-separated file where the 3rd element of each line
    is the nodes (separated by commas) expected for that validation range'''
    
    
if __name__ == "__main__":

    #loading parameters from command line that will be used below
    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    strainStr = sys.argv[2]
    chrom = sys.argv[3]
    validationFile = sys.argv[4]

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])

    strainNum = strainStr2NumD[strainStr]

    node = strainStr2NumD[strainStr]
    
    # load islands and genes
    islandByNodeL=islands.readIslands(paramD['islandOutFN'],tree,strainStr2NumD)

    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    
    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    familyT = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])

