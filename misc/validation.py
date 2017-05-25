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
    potentialIslands = islandsOnAllValidationNodes()
    returnIslands = []
    
    #for each island, check that it meets our criteria
    for island in potentialIslands:

        islandGenesInStrainL = getIslandGenesInStrain(island,strainNum,familyT)
        
        #check island's length > minLength
        if getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0) is not None:
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
    validationRanges, totalBases=readRanges()
    returnIslands = []  #holds islands that have any overlap
    overlapList=[] #holds the # of bp that overlap w/ a validation island for each island we check
    islandsList=[0]*len(validationRanges) #holds the number of xenoGI islands that overlap w/each validation range
    nodesLL,uniqueStrains = nodesPerRange()
    
    for island in potentialIslands:
        #get the start and end position for the islands
        islandGenesInStrainL = getIslandGenesInStrain(island,strainNum,familyT)
        if getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0) is not None:
            neighbGenesL,firstIslandGene,lastIslandGene=getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0)
            startPos = min(int(geneInfoD[geneNames.numToName(firstIslandGene)][4]), int(geneInfoD[geneNames.numToName(firstIslandGene)][5]))
            endPos = max(int(geneInfoD[geneNames.numToName(lastIslandGene)][5]),int(geneInfoD[geneNames.numToName(lastIslandGene)][4]))
            islandNode = island.mrca

            #check that island is in validation range
            inRange,overlap,indices=islandInRange(validationRanges, startPos, endPos, islandNode, nodesLL)

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

def islandInRange(validationRanges, startPos, endPos, islandNode, nodesLL):
    '''returns true if the island overlaps with any 
    validation ranges and false otherwise'''
    overlap = [0]*len(validationRanges)
    indices = [0]*len(validationRanges)
    for index in range(0,len(validationRanges)):
        vRange=validationRanges[index]
        if (islandNode in nodesLL[index]):
            overlap[index]=min((vRange[1]-vRange[0]),max((endPos-vRange[0]),0))-min((vRange[1]-vRange[0]),max((startPos-vRange[0]),0))
        if overlap[index]>0: indices[index]=1
    if sum(overlap)>0: return True, overlap, indices
    return False, overlap, indices


def islandsOnAllValidationNodes():
    '''makes a list of all islands at the nodes of interest'''
    nodesLL,uniqueStrains = nodesPerRange()
    allIslands = []
    for strain in uniqueStrains:
        allIslands+=islandByNodeL[int(strainStr2NumD[strain])]
    return allIslands

def readRanges():
    '''Read a file with start end ranges separated by tabs, and return list of lists.'''
    totalBases = 0
    f=open(validationFile,'r')
    outLL = []
    while True:
        s=f.readline()
        if s == "":
            break
        L=[]
        for elem in s.split('\t')[1:3]:
            L.append(int(elem))
        totalBases += (L[1]-L[0])
        outLL.append(L)
    return outLL, totalBases

def nodesPerRange():
    '''read a tab-separated file where the 3rd element of each line
    contains the strains (separated by commas) expected for that validation range'''
    f=open(validationFile,'r')
    nodesLL = []
    uniqueStrains = []
    while True:
        s=f.readline()
        if s == "":
            break
        nodes = []
        for strain in s.split('\t')[0].split(','):
            nodes.append(int(strainStr2NumD[strain]))
            if strain not in uniqueStrains: uniqueStrains.append(strain)
        nodesLL.append(nodes)
    return nodesLL,uniqueStrains

if __name__ == "__main__":

    #loading parameters from command line that will be used below
    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    strainStr = sys.argv[2]
    chrom = sys.argv[3]
    validationFile = sys.argv[4]

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])

    strainNum = strainStr2NumD[strainStr]

    #node = strainStr2NumD[strainStr]
    
    # load islands and genes
    islandByNodeL=islands.readIslands(paramD['islandOutFN'],tree,strainStr2NumD)

    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    
    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    familyT = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])

