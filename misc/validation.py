import sys,os,copy
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,genomes,trees,families,scores,islands,analysis

def islandsOfInterest():
    '''using input from the command line, prints out information about xenoGI's findings for the given validation ranges'''
    longIslands = islandsInStrainLongEnough(minGenes)
    longOnChrom = islandsOnChromosome(longIslands)
    overlapList,totalBases,islandsList,validationRanges,islandsPerRangeLL,coveragePerRangeL = islandsInRange(longOnChrom)
    overlap = sum(overlapList)
    for rangeIndex in range(0,len(validationRanges)):
        print(str(rangeIndex+1)+".","Range:",chrom,":",validationRanges[rangeIndex][0],"-",validationRanges[rangeIndex][1])
        covVal = (overlapList[rangeIndex])/(validationRanges[rangeIndex][1]-validationRanges[rangeIndex][0]) #calculate coverage for current range
        print("Coverage:",format(covVal,".3f"))
        print("Islands:",islandsPerRangeLL[rangeIndex])
        print("----")
    print("SUMMARY:")
    print("Ranges:",validationRanges)
    print("Number of Islands per range: ", islandsList)
    print("Percent overlap: ", format(overlap/totalBases,".3f"))
    print("All islands per range:",islandsPerRangeLL)

def islandsInStrainLongEnough(minGenes):
    '''returns a list of the xenoGI islands that are longer than minLength'''
    #list of islands in strain, empty list of islands to return
    potentialIslands = islandsOnAllValidationNodes()
    returnIslands = []
    
    #for each island, check that it meets our criteria
    for island in potentialIslands:
        islandGenesInStrainL = analysis.getIslandGenesInStrain(island,strainNum,familyL)
        
        #check island has at least min genes, if so add to list of potential islands
        if len(islandGenesInStrainL)>=minGenes: returnIslands.append(island)
    return returnIslands
    
def islandsOnChromosome(potentialIslands):
    '''of the potential islands passed in, returns a list of the islands on the desired 
    chromosome'''
    returnIslands = []
    #loop through each island,if an island is on the correct chromosome,
    #add it to our list of potential islands
    for island in potentialIslands:
        islandGenesInStrainL = analysis.getIslandGenesInStrain(island,strainNum,familyL)
        chromFound = geneInfoD[geneNames.numToName(islandGenesInStrainL[0])][3]
        if (chromFound == chrom): returnIslands.append(island)
    return returnIslands

def islandsInRange(potentialIslands):
    '''returns information about the islands in each validation range, including a list of islands in each range,
    coverage for each range, and values necessary to calculate total coverage'''
    validationRanges, totalBases=readRanges()
    islandsPerRangeLL = [[]]*len(validationRanges)
    #CoveragePerRangeL holds the start and end possitions of every island in each range
    coveragePerRangeLL = [[[],[]] for x in range(len(validationRanges))]
    #finalCoveragePerRangeL holds the total range covered by islands for each validation range
    finalCoveragePerRangeLL = [[0,0]]*len(validationRanges)
    returnIslands = []  #holds islands that have any overlap with any island
    overlapList=[] #holds the # of bp that overlap w/ a validation island for each island we check
    islandsList=[0]*len(validationRanges) #holds the number of xenoGI islands that overlap w/ each validation range
    nodesLL,uniqueStrains = nodesPerRange()
    
    for island in potentialIslands:
        #get the start and end position for the islands
        islandGenesInStrainL = analysis.getIslandGenesInStrain(island,strainNum,familyL)
        if analysis.getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0) is not None:
            neighbGenesL,firstIslandGene,lastIslandGene=analysis.getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0)
            startPos = min(int(geneInfoD[geneNames.numToName(firstIslandGene)][4]), int(geneInfoD[geneNames.numToName(firstIslandGene)][5]))
            endPos = max(int(geneInfoD[geneNames.numToName(lastIslandGene)][5]),int(geneInfoD[geneNames.numToName(lastIslandGene)][4]))
            islandNode = island.mrca

            #check that island is in validation range, if so, add it to the list of islands for that range
            # and add its start/end positions to coveragePerRangeLL
            inRange,indices=islandInRange(validationRanges, startPos, endPos, islandNode, nodesLL)
            for index in range(0,len(indices)):
                if indices[index] is 1:
                    islandsPerRangeLL[index]=islandsPerRangeLL[index]+[island.id]
                    coveragePerRangeLL[index][0].append(startPos)
                    coveragePerRangeLL[index][1].append(endPos)
                    

            #update the islandsList so it reflects how many xenoGI islands are in each range
            islandsList=list(map(lambda x,y:x+y, islandsList, indices))
               
    
            #add the island to the returnIslands list if it overlaps with any validation range
            if inRange:
                returnIslands.append(island)

    #update overlap list and finalCoveragePerRangeLL based on coveragePerRangeLL
    for valIndex in range(0,len(validationRanges)):
        finalCoveragePerRangeLL[valIndex][0] = min(coveragePerRangeLL[valIndex][0])
        if finalCoveragePerRangeLL[valIndex][0]<validationRanges[valIndex][0]: finalCoveragePerRangeLL[valIndex][0]=validationRanges[valIndex][0]
        finalCoveragePerRangeLL[valIndex][1] = max(coveragePerRangeLL[valIndex][1])
        if finalCoveragePerRangeLL[valIndex][1]>validationRanges[valIndex][1]: finalCoveragePerRangeLL[valIndex][1]=validationRanges[valIndex][1]

        overlapList.append(finalCoveragePerRangeLL[valIndex][1] - finalCoveragePerRangeLL[valIndex][0])
    return overlapList, totalBases, islandsList,validationRanges,islandsPerRangeLL,finalCoveragePerRangeLL

def islandInRange(validationRanges, startPos, endPos, islandNode, nodesLL):
    '''returns true if the island overlaps with any 
    validation ranges and false otherwise'''
    overlapLL = [[0,0]]*len(validationRanges)
    indices = [0]*len(validationRanges)
    for index in range(0,len(validationRanges)):
        vRange=validationRanges[index]
        if (islandNode in nodesLL[index]):
            overlapLL[index]=[min(vRange[1],max(vRange[0],startPos)), min(vRange[1],max(vRange[0],endPos))]
            if overlapLL[index][0] is vRange[1]:
                overlapLL[index][0]=vRange[0]
                overlapLL[index][1]=vRange[0]
        else: overlapLL[index]=[vRange[0],vRange[0]]
        if (overlapLL[index][1]-overlapLL[index][0])>0: indices[index]=1
    if sum(indices)>0: return True, overlapLL, indices
    return False, indices


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
    minGenes = int(sys.argv[5])

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])

    strainNum = strainStr2NumD[strainStr]

    #node = strainStr2NumD[strainStr]
    
    # load islands and genes
    islandByNodeL=islands.readIslands(paramD['islandOutFN'],tree,strainStr2NumD)

    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    
    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    familyL = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)

    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])

    islandsOfInterest()
