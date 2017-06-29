import sys,os,copy
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,genomes,trees,families,scores,islands,analysis

def islandsOfInterest():
    '''using input from the command line, prints out information about xenoGI's findings for the given validation ranges'''
    longIslands = islandsInStrainLongEnough(minGenes)
    longOnChromD = islandsOnChromosome(longIslands)
    overlapList,islandsList,islandsPerRangeLL,extra= islandsInRange(longOnChromD)
    overlap = sum(overlapList)
    for valRangeIndex in range(len(allValRanges)):
        print(str(valRangeIndex+1)+".","Range:",chromsL[valRangeIndex],":",allValRanges[valRangeIndex][0],"-",allValRanges[valRangeIndex][1])
        covVal = (overlapList[valRangeIndex])/(allValRanges[valRangeIndex][1]-allValRanges[valRangeIndex][0]) #calculate coverage for current range
        extraVal = (extra[valRangeIndex])/(allValRanges[valRangeIndex][1]-allValRanges[valRangeIndex][0]) #calculate percentage of extra coverage for current range
        print("Coverage:",format(covVal,".3f"))
        print("Extra Coverage:",format(extraVal,".3f"))
        print("Islands:",islandsPerRangeLL[valRangeIndex])
        print("----")
    print("SUMMARY:")
    print("Ranges:",allValRanges)
    print("Number of Islands per range: ", islandsList)
    print("Prop. ranges with 1 island: ",format(islandsList.count(1)/len(islandsList),".3f"))
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
    returnIslandsD = {}
    #loop through each island,if an island is on the correct chromosome,
    #add it to our list of potential islands
    for island in potentialIslands:
        islandGenesInStrainL = analysis.getIslandGenesInStrain(island,strainNum,familyL)
        chromFound = geneInfoD[geneNames.numToName(islandGenesInStrainL[0])][3]
        if chromFound in chromsL:
            if chromFound in returnIslandsD: returnIslandsD[chromFound].append(island)
            else: returnIslandsD[chromFound] = [island]
    return returnIslandsD

def islandsInRange(potentialIslandsD):
    '''returns information about the islands in each validation range, including a list of islands in each range,
    coverage for each range, and values necessary to calculate total coverage'''
    islandsPerRangeLL = [[]]*len(allValRanges)
    #CoveragePerRangeL holds the start and end possitions of every island in each range
    coveragePerRangeLL = [[[],[]] for x in range(len(allValRanges))]
    returnIslands = []  #holds islands that have any overlap with any island
    islandsList=[0]*len(allValRanges) #holds the number of xenoGI islands that overlap w/ each validation range

    for chrom,islandsL in potentialIslandsD.items():
        validationRangesForCurrChrom = []
        fillerRange = [0,0]
        for index in range(len(chromsL)):
            if chromsL[index]==chrom: validationRangesForCurrChrom.append(allValRanges[index])
            else:validationRangesForCurrChrom.append(fillerRange)
        for island in islandsL:
            #get the start and end position for the islands
            islandGenesInStrainL = analysis.getIslandGenesInStrain(island,strainNum,familyL)
            if analysis.getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0) is not None:
                neighbGenesL,firstIslandGene,lastIslandGene=analysis.getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,0)
                startPos = min(int(geneInfoD[geneNames.numToName(firstIslandGene)][4]), int(geneInfoD[geneNames.numToName(firstIslandGene)][5]))
                endPos = max(int(geneInfoD[geneNames.numToName(lastIslandGene)][5]),int(geneInfoD[geneNames.numToName(lastIslandGene)][4]))
                islandNode = island.mrca

                #check that island is in validation range, if so, add it to the list of islands for that range
                # and add its start/end positions to coveragePerRangeLL
                inRange,indices=islandInRange(validationRangesForCurrChrom, startPos, endPos, islandNode, nodesLL)
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
    overlapList,extra=overlapHelper(coveragePerRangeLL,allValRanges)
    return overlapList,islandsList,islandsPerRangeLL,extra

def overlapHelper(coveragePerRangeLL, valRanges):
    '''returns the number of bps covered for each validation range in a list and a the number of bps extra covered for each validation range'''
    overlap = [0]*len(valRanges)
    extra = [0]*len(valRanges)
    for valIndex in range(0,len(valRanges)):
        if len(coveragePerRangeLL[valIndex][0]) is not 0:
            rangeLength = valRanges[valIndex][1]-valRanges[valIndex][0]
            basesCoveredL = [0]*rangeLength
            for islandIndex in range(0,len(coveragePerRangeLL[valIndex][0])):
                #get the indices of the start and end of the island, if the island extends beyond the validation ranges,
                #reset that point to be the valRange start or end
                islandStart = max(coveragePerRangeLL[valIndex][0][islandIndex],valRanges[valIndex][0])-valRanges[valIndex][0]
                islandEnd = min(coveragePerRangeLL[valIndex][1][islandIndex],valRanges[valIndex][1])-valRanges[valIndex][0]
                for i in range(islandStart,islandEnd): basesCoveredL[i]=True
            overlap[valIndex]=basesCoveredL.count(True)
            extra[valIndex]=(valRanges[valIndex][0]-min(coveragePerRangeLL[valIndex][0]+[valRanges[valIndex][0]]))+(max(coveragePerRangeLL[valIndex][1]+[valRanges[valIndex][1]])-valRanges[valIndex][1])
    return overlap, extra
        

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
    if sum(indices)>0: return True,indices
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
    chromsL = []
    allValRanges = []
    f=open(validationFile,'r')
    while True:
        s=f.readline()
        if s == "":
            break

        chrom = s.split('\t')[0]
        chromsL.append(chrom)
        L=[]
        for elem in s.split('\t')[2:4]:
            L.append(int(elem))
        totalBases += (L[1]-L[0])
        allValRanges.append(L)
    return totalBases, chromsL, allValRanges

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
        for strain in s.split('\t')[1].split(','):
            nodes.append(int(strainStr2NumD[strain]))
            if strain not in uniqueStrains: uniqueStrains.append(strain)
        nodesLL.append(nodes)
    return nodesLL,uniqueStrains

if __name__ == "__main__":

    #loading parameters from command line that will be used below
    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    strainStr = sys.argv[2]
    validationFile = sys.argv[3]
    minGenes = int(sys.argv[4])

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

    totalBases, chromsL, allValRanges = readRanges()
    nodesLL,uniqueStrains = nodesPerRange()

    islandsOfInterest()
