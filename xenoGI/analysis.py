import sys
from .Family import *
from .Island import *
from . import trees
from . import scores
from . import islands

#### Analysis functions

## general

def printTable(L,indent=0,fileF=sys.stdout):
    '''Given tabular data in a list of lists (where sublists are rows)
print nicely so columns line up. Indent is an optional number of blank spaces to put in front of each row.'''
    # get max width for each column
    colMax=[]
    for col in range(len(L[0])):
        mx=0
        for row in L:
            if len(row[col]) > mx:
                mx = len(row[col])
        colMax.append(mx)
    
    # print
    for row in L:
        for col in range(len(row)):
            row[col]=row[col]+' ' * (colMax[col]-len(row[col]))
            
    for row in L:
        printStr = " "*indent + " | ".join(row)
        print(printStr.rstrip(),file=fileF)

def matchFamilyIsland(geneInfoD,geneNames,gene2FamD,fam2IslandD,searchStr):
    '''Return the island number, family number, and gene name(s)
associated with searchStr in geneInfoD. Searches for a match in all
fields of geneInfoD.'''
    # find matching gene names
    geneMatchL=[]
    for geneName in geneInfoD:
        valueT=geneInfoD[geneName]
        for value in (geneName,)+valueT:
            if type(value)==str:
                if searchStr in value:
                    geneMatchL.append(geneName)
                    break

    # get family numbers and island numbers
    outL=[]
    for geneName in geneMatchL:
        geneNum = geneNames.nameToNum(geneName)
        fam=gene2FamD[geneNum]
        isl=fam2IslandD[fam]
        outL.append((geneName,fam,isl.id))
    return outL
        
## Print scores associated with a family

def printScoreMatrix(familyNum,subtreeL,familyL,geneNames,scoresO,scoreType,fileF):
    '''Print a matrix of scores between all the genes in a familyNum. Scores
are provided by scoresO, and we're extracting the values associated
with scoreType in the edges of this graph.
    '''

    familyGeneNumsL = familyL[familyNum].getGeneNums()
    
    rowsL = []
    geneNamesL = [geneNames.numToName(gn) for gn in familyGeneNumsL]
    rowsL.append([''] + geneNamesL)
    
    for rowi,gn1 in enumerate(familyGeneNumsL):
        row = [geneNames.numToName(familyGeneNumsL[rowi])]
        for gn2 in familyGeneNumsL:
            if scoresO.isEdgePresentByEndNodes(gn1,gn2):
                row.append(format(scoresO.getScoreByEndNodes(gn1,gn2,scoreType),".3f"))
            else:
                row.append('-')
        rowsL.append(row)

    printTable(rowsL,indent=2,fileF=fileF)

def printOutsideFamilyScores(familyNum,subtreeL,familyL,geneNames,scoresO,fileF):
    '''Given a family, print scores for all non-family members with a
connection to genes in family. Scores are provided in the network
scoresO.
    '''

    family = familyL[familyNum]
    outsideGeneNumsT = family.getOutsideConnections(scoresO)
    
    rowL = []
    for familyGeneNum in family.getGeneNums():
        familyGeneName = geneNames.numToName(familyGeneNum)
        for outsideGeneNum in outsideGeneNumsT:
            if scoresO.isEdgePresentByEndNodes(familyGeneNum,outsideGeneNum):
                outsideGeneName = geneNames.numToName(outsideGeneNum)
                rawSc=scoresO.getScoreByEndNodes(familyGeneNum,outsideGeneNum,'rawSc')
                normSc=scoresO.getScoreByEndNodes(familyGeneNum,outsideGeneNum,'normSc')
                coreSynSc=scoresO.getScoreByEndNodes(familyGeneNum,outsideGeneNum,'coreSynSc')
                synSc=scoresO.getScoreByEndNodes(familyGeneNum,outsideGeneNum,'synSc')
                rowL.append([familyGeneName,outsideGeneName,format(rawSc,".3f"),format(normSc,".3f"),format(coreSynSc,".3f"),format(synSc,".3f")])

    rowL.sort(key=lambda x: x[2],reverse=True) # sort by score
    rowL.insert(0,['----------','-----------','---','----','-------','---'])
    rowL.insert(0,['Inside fam','Outside fam','Raw','Norm','CoreSyn','Syn'])
                
    print("Printing all scores with non-family members",file=fileF)
    printTable(rowL,indent=2,fileF=fileF)

## Print strains where family is present, and not present

def familyPrintStrainsPresentAbsent(tree,strainNum2StrD,familyL,famNum,fileF=sys.stdout):
    '''Print a list of strains where the family is present, and another where it is absent.'''

    presL=[]
    notPresL=[]
    for leafNum in trees.leafList(tree):
        if familyL[famNum].isInStrain(leafNum):
            presL.append(strainNum2StrD[leafNum])
        else:
            notPresL.append(strainNum2StrD[leafNum])
    print("Family:",famNum,file=fileF)
    print("  Strains possessing:",file=fileF)
    for strain in presL:
        print("    "+strain,file=fileF)
    print(file=fileF)
    print("  Strains lacking:",file=fileF)
    for strain in notPresL:
        print("    "+strain,file=fileF)


    
## Print all islands at node

def printIslandLSummary(island,fileF):
    '''Given a list of islands in island (ie a list from a single node),
print a simple tabular summary indicating how many families they
have.
    '''
    lenL = []
    for isl in island:
        lenL.append(len(isl)) # len of island is num families

    # count how many times each length occurs
    lnCtD = {}
    for ln in lenL:
        if ln in lnCtD:
            lnCtD[ln] += 1
        else:
            lnCtD[ln] = 1

    # print out
    printL = []
    row = ['Num families in island','Number of occurrences']
    printL.append(row)
    
    for ln,occurrences in sorted(lnCtD.items()):
        printL.append([str(ln), str(occurrences)])

    printTable(printL,indent=8,fileF=fileF)
    
def vPrintIsland(island,subtreeL,familyL,strainNum2StrD,geneNames,geneInfoD,fileF):
    '''Verbose print of an island.'''

    print("  Island",island.id,file=fileF)
    
    # get species nodes subtended by this mrca
    speciesNodesL=trees.leafList(subtreeL[island.mrca])

    # put everything in lists.
    printL=[]
    printL.append(['Family'])
    for node in speciesNodesL:
        printL[0].append(strainNum2StrD[node])
    for fam in island.familyL:
        newRow=[]
        newRow.append(str(fam))
        for node in speciesNodesL:
            geneT = familyL[fam].famGeneT[node]
            for geneNum in geneT:
                geneName = geneNames.numToName(geneNum)
                commonGeneName = "("+geneInfoD[geneName][0]+")" if geneInfoD[geneName][0] != '' else ''
                newRow.append(",".join([geneName + commonGeneName for geneNum in geneT]))
        printL.append(newRow)
    printTable(printL,indent=4,fileF=fileF)


def vPrintIslands(islandL,subtreeL,familyL,strainNum2StrD,geneNames,geneInfoD,fileF):
    '''Print a list of islands.'''
    print("  Summary",file=fileF)
    printIslandLSummary(islandL,fileF)
    print("  ---",file=fileF)
    for island in islandL:
        vPrintIsland(island,subtreeL,familyL,strainNum2StrD,geneNames,geneInfoD,fileF)
        print('  ---',file=fileF)

def vPrintAllIslands(islandByNodeL,tree,rootFocalClade,subtreeL,familyL,strainStr2NumD,strainNum2StrD,geneNames,geneInfoD,fileF):
    '''Loop over all nodes in tree, printing islands at each. '''
    rootFocalCladeNum = strainStr2NumD[rootFocalClade]
    focalTree = trees.subtree(tree,rootFocalCladeNum)
    for node in trees.nodeList(focalTree):
        nodeStr = strainNum2StrD[node]
        print('########################### ',"Islands at node",nodeStr,file=fileF)
        print('',file=fileF)
        vPrintIslands(islandByNodeL[node],subtreeL,familyL,strainNum2StrD,geneNames,geneInfoD,fileF)

        
def createGene2FamD(familyL):
    '''Given the family information in familyL, create a dictionary
gene2FamD which maps from gene number to family number.'''
    gene2FamD={}
    for famNum in range(len(familyL)):
        for geneT in familyL[famNum].famGeneT:
            for gene in geneT:
                gene2FamD[gene]=famNum
    return gene2FamD

def createFam2IslandD(islandL):
    '''Given islandL, our list of islands, create a dictionary that maps
family number to island number.
    '''
    fam2IslandD={}
    for islandsAtNodeL in islandL:
        for island in islandsAtNodeL:
            for famNum in island.familyL:
                fam2IslandD[famNum]=island
    return fam2IslandD

## Print species files with all the genes, grouped by contig

def printSpeciesContigs(geneOrderT,fileStemStr,fileExtensionStr,geneNames,gene2FamD,fam2IslandD,geneInfoD,familyL,strainNum2StrD):
    '''This function produces a set of species specific genome
files. These contain all the genes in a strain laid out in the order
they occur on the contigs. Each gene entry include island and family
information, as well as a brief description of the gene's function.'''

    for strainNum in range(len(geneOrderT)):
        strainStr = strainNum2StrD[strainNum]
        contigT = geneOrderT[strainNum]
        if contigT != None:
            # it's a tip
            fileF=open(fileStemStr+'-'+strainStr+fileExtensionStr,'w')
            for contig in contigT:
                print("########### Contig",file=fileF)
                printGenes(contig,geneNames,gene2FamD,fam2IslandD,geneInfoD,[],familyL,strainNum2StrD,fileF)

            fileF.close()
    
## Print neighborhood of an island

def printIslandNeighb(islandNum,synWSize,subtreeL,islandByNodeL,familyL,geneOrderT,gene2FamD,fam2IslandD,geneInfoD,geneNames,strainNum2StrD,fileF):
    '''Print the neighborhood of an island. We include the genes in the island and synWSize/2 genes in either direction.'''

    print("  Island:",islandNum,file=fileF)
    
    genesInEitherDirec = int(synWSize/2)

    # get the island object for this islandNum
    for listOfIslands in islandByNodeL:
        _,island = islands.searchIslandsByID(listOfIslands,islandNum)
        if island != None: break

    if island == None:
        raise ValueError("Island "+str(islandNum)+" not found.")
        
    mrca = island.mrca
    print("  mrca:",strainNum2StrD[mrca],file=fileF)

    leavesL=trees.leafList(subtreeL[mrca])

    for strainNum in leavesL:

        print("  In",strainNum2StrD[strainNum],end=' ',file=fileF)

        islandGenesInStrainL = getIslandGenesInStrain(island,strainNum,familyL)

        if islandGenesInStrainL == []:
            print("the island is not found.",file=fileF)
        else:

            neighbGenesL,firstIslandGene,lastIslandGene=getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,genesInEitherDirec)

            # print coordinates of island in this strain
            chrom=geneInfoD[geneNames.numToName(islandGenesInStrainL[0])][3]
            startPos = geneInfoD[geneNames.numToName(firstIslandGene)][4]
            endPos = geneInfoD[geneNames.numToName(lastIslandGene)][5]

            print("(Coordinates",chrom+":"+str(startPos)+"-"+str(endPos)+")",file=fileF)

            printGenes(neighbGenesL,geneNames,gene2FamD,fam2IslandD,geneInfoD,islandGenesInStrainL,familyL,strainNum2StrD,fileF)


def getIslandGenesInStrain(island,strainNum,familyL):
    '''Given an island, a strain number, and our tuple of family
objects, return all the genes in the island for that strain.'''
    genesL=[]
    for familyNum in island.familyL:
        geneT=familyL[familyNum].famGeneT[strainNum]
        genesL.extend(geneT)
    return genesL

def getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,genesInEitherDirec):
    ''''''
    neighbGenesL=[]
    for contig in geneOrderT[strainNum]:
        try:
            # get index of all of these. We're assuming they're on
            # the same contig and in the same area.
            indL=[contig.index(gene) for gene in islandGenesInStrainL]
            maxInd = max(indL)
            minInd = min(indL)
            end = maxInd + genesInEitherDirec +1
            st = minInd-genesInEitherDirec if minInd-genesInEitherDirec>0 else 0 # st can't be less than 0
            neighbGenesL=contig[st:end]

            # get gene numbers of first and last genes in island
            firstGene = contig[minInd]
            lastGene = contig[maxInd]
            
            return neighbGenesL,firstGene,lastGene
        except ValueError:
            continue

def printGenes(neighbGenesL,geneNames,gene2FamD,fam2IslandD,geneInfoD,islandGenesInStrainL,familyL,strainNum2StrD,fileF):
    '''Given a list of contiguous genes, print them out nicely with information on family and island etc. Put * next to any that are in islandGenesInStrainL.'''
            
    # now print the neighbors
    rowsL=[]
    for tempGene in neighbGenesL:
        tempGeneName=geneNames.numToName(tempGene)
        tempFamNum=gene2FamD[tempGene]
        tempGeneIsland=fam2IslandD[tempFamNum]

        if tempGeneName in geneInfoD:
            descrip = geneInfoD[tempGeneName][2]
        else:
            descrip = ''

        # mark genes in the island with a *
        if tempGene in islandGenesInStrainL:
            tempGeneName = '* '+tempGeneName
        else:
            tempGeneName = '  '+tempGeneName

        infoL = [tempGeneName,"isl:"+str(tempGeneIsland.id),"fam:"+str(tempFamNum),"errSc:"+str(familyL[tempFamNum].possibleErrorCt),"mrca:"+strainNum2StrD[tempGeneIsland.mrca],descrip]

        rowsL.append(infoL)

    printTable(rowsL,indent=4,fileF=fileF)

# support functions for printCoreNonCoreByNode

def createFamilyByNodeL(geneOrderT,gene2FamD):
    '''Create a list where index is node, and where the value is all the
families present in that strain. Internal nodes are None, tips are a set of
families.'''

    # an index for each node. internal nodes get none, just like in geneOrderT
    familyByNodeL=[set() if x!=None else None for x in geneOrderT]

    for i in range(len(geneOrderT)):
        contigT = geneOrderT[i]
        if contigT != None:
            for contig in contigT:
                for geneNum in contig:
                    fam=gene2FamD[geneNum]
                    familyByNodeL[i].add(fam)
    return familyByNodeL

def coreNonCoreCtAtNode(tree,node,familyByNodeL,familyL):
    '''Given a tree and a node, first get all the families present in
descendant species. Then figure out which of these families are
non-core (their mrca is located below node) and which are core (mrca
is at node or above). Return count of non-core and core.'''

    subtree = trees.subtree(tree,node)
    nonCoreMrcaL = trees.nodeList(subtree[1]) + trees.nodeList(subtree[2])
    
    # get set of all families with members in descendant species of this node
    decFamS = set()
    for leaf in trees.leafList(subtree):
        decFamS.update(familyByNodeL[leaf])

    # figure out which are core, non-core
    coreCt=0
    nonCoreCt=0
    totCt = len(decFamS)
    for fam in decFamS:
        if familyL[fam].mrca in nonCoreMrcaL:
            nonCoreCt += 1
        else:
            coreCt += 1
            
    return nonCoreCt,coreCt
