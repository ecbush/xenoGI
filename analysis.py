from Family import *
from Island import *
import trees,scores,islands
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages

#### Analysis functions

## general

def printTable(L,indent=0):
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
        print(printStr.rstrip())


## Print scores associated with a family

def printScoreMatrix(family,subtreeL,familyT,geneNames,G,scoreType):
    '''Print a matrix of scores between all the genes in a family. Scores
are provided by the graph G, and we're extracting the values associated with scoreType in the edges of this graph.'''

    familyGenesL = getGenesInFamily(family,subtreeL,familyT)

    rowsL = []
    geneNamesL = [geneNames.numToName(gn) for gn in familyGenesL]
    rowsL.append([''] + geneNamesL)
    
    for rowi,gn1 in enumerate(familyGenesL):
        row = [geneNames.numToName(familyGenesL[rowi])]
        for gn2 in familyGenesL:
            data=G.get_edge_data(gn1,gn2)
            if data == None:
                row.append('-')
            else:
                row.append(format(data[scoreType],".3f"))
        rowsL.append(row)

    printTable(rowsL,2)

def printOutsideFamilyScores(family,subtreeL,familyT,geneNames,scoresG):
    '''Given a family, print scores for all non-family members with a
connection to genes in family. Scores are provided in the network
rawScoresG.
    '''
    
    familyGenesL = getGenesInFamily(family,subtreeL,familyT)

    rowL = []
    for gene in familyGenesL:
        geneName = geneNames.numToName(gene)
        for edge in scoresG.edges_iter(gene):
            if not edge[1] in familyGenesL:
                # this connection is with a gene outside the family
                otherGeneName = geneNames.numToName(edge[1])
                rawSc=scoresG.get_edge_data(gene,edge[1])['rawSc']
                normSc=scoresG.get_edge_data(gene,edge[1])['normSc']
                coreSynSc=scoresG.get_edge_data(gene,edge[1])['coreSynSc']
                synSc=scoresG.get_edge_data(gene,edge[1])['synSc']
                rowL.append([geneName,otherGeneName,format(rawSc,".3f"),format(normSc,".3f"),format(coreSynSc,".3f"),format(synSc,".3f")])

    rowL.sort(key=lambda x: x[2],reverse=True) # sort by score
    rowL.insert(0,['----------','-----------','---','----','-------','---'])
    rowL.insert(0,['Inside fam','Outside fam','Raw','Norm','CoreSyn','Syn'])
                
    print("Printing all scores with non-family members")
    printTable(rowL,2)




## Print all islands at node

def printIslandLSummary(island):
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

    printTable(printL,8)
    
def vPrintIsland(island,subtreeL,familyT,strainNum2StrD,geneNames):
    '''Verbose print of an island.'''

    print("  Island",island.id)
    
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
            ct,genesL = familyT[fam].famGeneT[node]
            newRow.append(",".join([geneNames.numToName(gene) for gene in genesL]))
        printL.append(newRow)
    printTable(printL,4)


def vPrintIslands(islandL,subtreeL,familyT,strainNum2StrD,geneNames):
    '''Print a list of islands.'''
    print("Summary of islands")
    printIslandLSummary(islandL)
    print("Print outs of each island")
    for island in islandL:
        vPrintIsland(island,subtreeL,familyT,strainNum2StrD,geneNames)
        print('  ---')

def createGene2FamD(familyT):
    '''Given the family information in familyT, create a dictionary
gene2FamD which maps from gene number to family number.'''
    gene2FamD={}
    for famNum in range(len(familyT)):
        for gnCt,geneT in familyT[famNum].famGeneT:
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

def getGenesInFamily(familyNum,subtreeL,familyT):
    '''Given a family, return a list of all genes in it (in numerical form).'''
    mrca = familyT[familyNum].mrca
    leavesL=trees.leafList(subtreeL[mrca])

    genesL=[]
    for strainNum in leavesL:

        geneT=familyT[familyNum].famGeneT[strainNum][1]
        genesL.extend(geneT)

    return genesL


## Print neighborhood of an island

def getIslandGenesinStrain(island,strainNum,familyT):
    '''Given an island, a strain number, and our tuple of family
objects, return all the genes in the island for that strain.'''
    genesL=[]
    for familyNum in island.familyL:
        geneT=familyT[familyNum].famGeneT[strainNum][1]
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
            end = max(indL) + genesInEitherDirec +1
            st = min(indL)-genesInEitherDirec if min(indL)-genesInEitherDirec>0 else 0 # st can't be less than 0
            neighbGenesL=contig[st:end]
            return neighbGenesL
        except ValueError:
            continue

    
def printIslandNeighb(islandNum,synWSize,subtreeL,islandByNodeL,familyT,geneOrderT,gene2FamD,fam2IslandD,geneInfoD,geneNames,strainNum2StrD):
    '''Print the neighborhood of an island. We include the genes in the island and synWSize/2 genes in either direction.'''

    genesInEitherDirec = int(synWSize/2)

    # get the island object for this islandNum
    for listOfIslands in islandByNodeL:
        _,island = islands.searchIslandsByID(listOfIslands,islandNum)
        if island != None: break

    mrca = island.mrca
    print("  Island mrca:",strainNum2StrD[mrca])

    leavesL=trees.leafList(subtreeL[mrca])
    print("  Island found in the following species:")
    for strainNum in leavesL:
        print("    ",strainNum2StrD[strainNum])
    
    for strainNum in leavesL:

        print("  Neighbors for island",islandNum,"in",strainNum2StrD[strainNum])

        islandGenesInStrainL = getIslandGenesinStrain(island,strainNum,familyT)
        
        if islandGenesInStrainL == []:
            print("    There are no island members in this species.")

        neighbGenesL=getNeighborhoodGenes(strainNum,geneOrderT,islandGenesInStrainL,genesInEitherDirec)
            
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
                
            infoL = [tempGeneName,"isl:"+str(tempGeneIsland.id),"fam:"+str(tempFamNum),"mrca:"+strainNum2StrD[tempGeneIsland.mrca],descrip]

            rowsL.append(infoL)

        printTable(rowsL,4)

    
def printFamNeighb(family,synWSize,subtreeL,familyT,geneOrderT,gene2FamD,fam2IslandD,geneInfoD,geneNames,strainNum2StrD):
    '''Family Neighborhood. Given a family, print out
the families nearby in each of its species and going either direction,
within synWSize/2.
    '''

    genesInEitherDirec = int(synWSize/2)
    
    mrca = familyT[family].mrca
    print("  Family mrca:",strainNum2StrD[mrca])

    leavesL=trees.leafList(subtreeL[mrca])
    print("  Tree includes")
    for leaf in leavesL:
        print("    ",strainNum2StrD[leaf])
    
    for leaf in leavesL:

        print("  Neighbors for family",family,"in",strainNum2StrD[leaf])

        geneT=familyT[family].famGeneT[leaf][1]
        if geneT == ():
            print("    There are no family members in this species.")

        for gene in geneT:

            print("    For",geneNames.numToName(gene))
            
            neighbGenesL=[]
            for contig in geneOrderT[leaf]:
                try:
                    ind=contig.index(gene)
                    end = ind + genesInEitherDirec
                    st = ind-genesInEitherDirec if ind-genesInEitherDirec>0 else 0 # st can't be less than 0
                    neighbGenesL=contig[st:end]
                    break
                except ValueError:
                    continue

            rowsL=[]
            for tempGene in neighbGenesL:
                geneName=geneNames.numToName(tempGene)
                famNum=gene2FamD[tempGene]
                island=fam2IslandD[famNum]

                if geneName in geneInfoD:
                    descrip = geneInfoD[geneName][2]
                else:
                    descrip = ''

                # if its the one in the family we're querying, mark with *
                if tempGene == gene:
                    geneName = '*'+geneName
                else:
                    geneName = ' '+geneName
                infoL = [geneName,"fam:"+str(famNum),"gr:"+str(island.id),"mrca:"+strainNum2StrD[island.mrca],descrip]

                rowsL.append(infoL)

            printTable(rowsL,4)



# Plots of scores

def scoreHists(scoresFN,outFN,numBins,geneNames,scoreType):
    '''Read through a scores file, and separate into all pairwise comparisons. Then plot hist of each.'''

    # currently, this seems to require a display for interactive
    # plots. would be nice to make it run without that...

    pairD = readScorePairs(scoresFN,geneNames,scoreType)

    pyplot.ioff() # turn off interactive mode
    with PdfPages(outFN) as pdf:
        for key in pairD:
            fig = pyplot.figure()
            pyplot.hist(pairD[key],bins=numBins)
            pyplot.title('-'.join(key))
            pdf.savefig()
            pyplot.close()

    #pyplot.show()
            

def readScorePairs(scoresFN,geneNames,scoreType):
    '''Read through a scores file, and separate into all pairwise
comparisons. Return as dict.'''
    
    pairD = {}

    scoreG = scores.readGraph(scoresFN,geneNames=None)
    
    for gn1,gn2,sc in scoreG.edges_iter(data=scoreType):

        sp1 = geneNames.numToStrainName(gn1)
        sp2 = geneNames.numToStrainName(gn2)
        key = tuple(sorted([sp1,sp2]))
        
        if key in pairD:
            pairD[key].append(sc)
        else:
            pairD[key] = [sc]
        
    return pairD
