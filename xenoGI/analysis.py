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

def matchFamilyIsland(genesO,gene2FamIslandD,searchStr):
    '''Return the island number, family number, and gene name(s)
associated with searchStr in genesO.geneInfoD. Searches for a match in all
fields of geneInfoD.'''
    # find matching gene names
    geneMatchL=[]
    for geneNum in genesO.geneInfoD:
        valueT=genesO.numToGeneInfo(geneNum)
        for value in valueT:
            if type(value)==str:
                if searchStr in value:
                    geneMatchL.append(geneNum)
                    break

    # get family numbers and island numbers
    outL=[]
    for geneNum in geneMatchL:
        geneName = genesO.numToName(geneNum)
        (locusIslandNum, famNum, locusFamNum) = gene2FamIslandD[geneNum]
        outL.append((geneName,locusIslandNum, famNum, locusFamNum))
    return outL
        
## Print scores associated with a family

def printScoreMatrix(familyNum,familiesO,genesO,scoresO,scoreType,fileF):
    '''Print a matrix of scores between all the genes in a family given by
familyNum. Scores are provided by scoresO, and we're extracting the
values associated with scoreType in the edges of this graph.
    '''

    familyGeneNumsL = []
    for lfO in familiesO.getFamily(familyNum).getLocusFamilies():
        familyGeneNumsL.extend(lfO.iterGenes())
    
    rowsL = []
    geneNamesL = [genesO.numToName(gn) for gn in familyGeneNumsL]
    rowsL.append([''] + geneNamesL)
    
    for rowi,gn1 in enumerate(familyGeneNumsL):
        row = [genesO.numToName(familyGeneNumsL[rowi])]
        for gn2 in familyGeneNumsL:
            if scoresO.isEdgePresentByEndNodes(gn1,gn2):
                row.append(format(scoresO.getScoreByEndNodes(gn1,gn2,scoreType),".3f"))
            else:
                row.append('-')
        rowsL.append(row)

    printTable(rowsL,indent=2,fileF=fileF)

def printOutsideFamilyScores(familyNum,familiesO,genesO,scoresO,fileF):
    '''Given a family, print scores for all non-family members with a
connection to genes in family. Scores are provided in the network
scoresO.
    '''

    family = familiesO.getFamily(familyNum)
    outsideGeneNumsS = family.getOutsideConnections(scoresO)
    
    rowL = []
    for familyGeneNum in family.iterGenes():
        familyGeneName = genesO.numToName(familyGeneNum)
        for outsideGeneNum in outsideGeneNumsS:
            if scoresO.isEdgePresentByEndNodes(familyGeneNum,outsideGeneNum):
                outsideGeneName = genesO.numToName(outsideGeneNum)
                rawSc=scoresO.getScoreByEndNodes(familyGeneNum,outsideGeneNum,'rawSc')
                synSc=scoresO.getScoreByEndNodes(familyGeneNum,outsideGeneNum,'synSc')
                coreSynSc=scoresO.getScoreByEndNodes(familyGeneNum,outsideGeneNum,'coreSynSc')
                rowL.append([familyGeneName,outsideGeneName,format(rawSc,".3f"),format(synSc,".3f"),format(coreSynSc,".3f")])

    rowL.sort(key=lambda x: x[2],reverse=True) # sort by score
    rowL.insert(0,['----------','-----------','---','---','-------'])
    rowL.insert(0,['Inside fam','Outside fam','Raw','Syn','CoreSyn'])
                
    print("Printing all scores with non-family members",file=fileF)
    printTable(rowL,indent=2,fileF=fileF)

    
## Print all islands at node

def printIslandLSummary(islandL,fileF):
    '''Given a list of locus islands in islandL (ie a list from a single node),
print a simple tabular summary indicating how many locus families they
have.
    '''
    lenL = []
    for isl in islandL:
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
    row = ['Num LocusFamilies in LocusIsland','Number of occurrences']
    printL.append(row)
    
    for ln,occurrences in sorted(lnCtD.items()):
        printL.append([str(ln), str(occurrences)])

    printTable(printL,indent=8,fileF=fileF)

def vPrintLocusIsland(island,subtreeD,familiesO,genesO,fileF):
    '''Verbose print of a locus island.'''

    print("  LocusIsland",island.id,file=fileF)
    
    # get species nodes subtended by this mrca
    speciesNodesL=trees.leafList(subtreeD[island.mrca])

    # put everything in lists.
    printL=[]
    printL.append(['LocusFamily','Family'])
    for node in speciesNodesL:
        printL[0].append(node)
        
    for locusFamO in island.iterLocusFamilies(familiesO):
        newRow=[]
        newRow.append(str(locusFamO.locusFamNum))
        newRow.append(str(locusFamO.famNum))
        for node in speciesNodesL:
            entryL = []
            for geneNum in locusFamO.iterGenesByStrain(node):
                infoT=genesO.numToGeneInfo(geneNum)
                geneName,commonName,locusTag,descrip,chrom,start,end,strand=infoT
                commonGeneName = "("+commonName+")" if commonName != '' else ''
                entryL.append(geneName + commonGeneName)
                
            if entryL == []:
                entry = ''
            else:
                entry = ",".join(entryL)

            newRow.append(entry)
        printL.append(newRow)
       
    printTable(printL,indent=4,fileF=fileF)

def vPrintLocusIslandsAtNode(islandL,subtreeD,familiesO,genesO,fileF):
    '''Print a list of islands at a single node.'''
    print("  Summary",file=fileF)
    printIslandLSummary(islandL,fileF)
    print("  ---",file=fileF)
    for island in islandL:
        vPrintLocusIsland(island,subtreeD,familiesO,genesO,fileF)
        print('  ---',file=fileF)

def vPrintAllLocusIslands(islandByNodeD,tree,rootFocalClade,subtreeD,familiesO,genesO,fileF):
    '''Loop over all nodes in tree, printing islands at each. '''
    focalTree = trees.subtree(tree,rootFocalClade)
    for node in trees.nodeList(focalTree):
        print('########################### ',"Locus Islands at node",node,file=fileF)
        print('',file=fileF)
        vPrintLocusIslandsAtNode(islandByNodeD[node],subtreeD,familiesO,genesO,fileF)

def createGene2FamIslandD(islandByNodeD,familiesO):
    '''Creates a dictionary keyed by gene number which has the
LocusIsland, Family and LocusFamily for each gene.'''
    D = {}
    for islandsAtNodeL in islandByNodeD.values():
        for locusIslandO in islandsAtNodeL:
            for locusFamO in locusIslandO.iterLocusFamilies(familiesO):
                for gene in locusFamO.iterGenes():
                    D[gene] = (locusIslandO.id,locusFamO.famNum,locusFamO.locusFamNum)

    return D

## Print species files with all the genes, grouped by contig

def printSpeciesContigs(geneOrderD,fileStemStr,fileExtensionStr,genesO,gene2FamIslandD,familiesO,strainNamesT):
    '''This function produces a set of species specific genome
files. These contain all the genes in a strain laid out in the order
they occur on the contigs. Each gene entry includes LocusIsland and LocusFamily
information, as well as a brief description of the gene's function.'''

    for strainName in strainNamesT:
        if strainName in geneOrderD:
            contigT = geneOrderD[strainName]
            with open(fileStemStr+'-'+strainName+fileExtensionStr,'w') as fileF:
                for contig in contigT:
                    print("########### Contig",file=fileF)
                    printGenes(contig,genesO,gene2FamIslandD,[],familiesO,fileF)
    return
    
## Print neighborhood of an island

def printLocusIslandNeighb(islandNum,synWSize,subtreeD,islandByNodeD,familiesO,geneOrderD,gene2FamIslandD,genesO,fileF):
    '''Print the neighborhood of an island. We include the genes in the island and synWSize/2 genes in either direction.'''

    print("  LocusIsland:",islandNum,file=fileF)
    
    genesInEitherDirec = int(synWSize/2)

    # get the island object for this islandNum
    for listOfIslands in islandByNodeD.values():
        _,island = islands.searchLocIslandsByID(listOfIslands,islandNum)
        if island != None: break

    if island == None:
        raise ValueError("LocusIsland "+str(islandNum)+" not found.")
        
    mrca = island.mrca
    print("  mrca:",mrca,file=fileF)

    leavesL=trees.leafList(subtreeD[mrca])

    for strainName in leavesL:

        print("  In",strainName,end=' ',file=fileF)

        islandGenesInStrainL = getIslandGenesInStrain(island,strainName,familiesO)

        if islandGenesInStrainL == []:
            print("the island is not found.",file=fileF)
        else:

            neighbGenesL,firstIslandGene,lastIslandGene=getNeighborhoodGenes(strainName,geneOrderD,islandGenesInStrainL,genesInEitherDirec)

            # print coordinates of island in this strain
            chrom = genesO.numToGeneInfo(islandGenesInStrainL[0])[4]
            startPos = genesO.numToGeneInfo(firstIslandGene)[5]
            endPos = genesO.numToGeneInfo(lastIslandGene)[6]
            
            print("(Coordinates",chrom+":"+str(startPos)+"-"+str(endPos)+")",file=fileF)
            printGenes(neighbGenesL,genesO,gene2FamIslandD,islandGenesInStrainL,familiesO,fileF)


def getIslandGenesInStrain(island,strainName,familiesO):
    '''Given an island, a strain number, and our families object, return
all the genes in the island for that strain.'''
    genesL=[]
    for locFam in island.iterLocusFamilies(familiesO):
        genesL.extend(locFam.iterGenesByStrain(strainName))
    return genesL

def getNeighborhoodGenes(strainName,geneOrderD,islandGenesInStrainL,genesInEitherDirec):
    ''''''
    neighbGenesL=[]
    for contig in geneOrderD[strainName]:
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

def printGenes(neighbGenesL,genesO,gene2FamIslandD,islandGenesInStrainL,familiesO,fileF):
    '''Given a list of contiguous genes, print them out nicely with
information on locus family and locus island etc. Put * next to any
that are in islandGenesInStrainL.'''

    # now print the neighbors
    rowsL=[]
    for geneNum in neighbGenesL:

        infoT=genesO.numToGeneInfo(geneNum)
        geneName,commonName,locusTag,descrip,chrom,start,end,strand=infoT

        
        locIslNum,famNum,locFamNum = gene2FamIslandD[geneNum]
        lfMrca = familiesO.getLocusFamily(locFamNum).lfMrca

        
        #if geneName in geneInfoD:
        #    descrip = geneInfoD[geneName][2]
        #else:
        #    descrip = ''

        # mark genes in the island with a *
        if geneNum in islandGenesInStrainL:
            geneName = '* '+geneName
        else:
            geneName = '  '+geneName

        infoL = [geneName,"locIsl:"+str(locIslNum),"fam:"+str(famNum),"locFam:"+str(locFamNum),"locFamMRCA:"+lfMrca,descrip]

        rowsL.append(infoL)

    printTable(rowsL,indent=4,fileF=fileF)

# support functions for printCoreNonCoreByNode

def createFamilyByNodeL(geneOrderT,gene2LocFamD):
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
                    fam=gene2LocFamD[geneNum]
                    familyByNodeL[i].add(fam)
    return familyByNodeL

def coreNonCoreCtAtNode(tree,node,familyByNodeL,familiesO):
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
