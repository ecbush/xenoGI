import sys
import genomes,trees,families,scores
from Family import *
from Group import *

## Load groups

def readGroupOut(fn,tree,strainStr2NumD):
    '''Given a file name for a groups output file, load back
recreating groupByNodeL.'''

    groupByNodeL=[[] for i in range(trees.nodeCount(tree))]
    
    f=open(fn,'r')
    while True:
        s=f.readline()
        if s == '':
            break
        gr=str2Group(s.rstrip(),strainStr2NumD)
        groupByNodeL[gr.mrca].append(gr)
    f.close()

    return groupByNodeL

## Analysis functions

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

def printGroupLSummary(group):
    '''Given a list of groups in group (ie a list from a single node),
print a simple tabular summary indicating how many families they
have.
    '''
    lenL = []
    for gr in group:
        lenL.append(len(gr)) # len of group is num families

    # count how many times each length occurs
    lnCtD = {}
    for ln in lenL:
        if ln in lnCtD:
            lnCtD[ln] += 1
        else:
            lnCtD[ln] = 1

    # print out
    printL = []
    row = ['Num families in group','Number of occurrences']
    printL.append(row)
    
    for ln,occurrences in sorted(lnCtD.items()):
        printL.append([str(ln), str(occurrences)])

    printTable(printL,8)
            
    
    
def vPrintGroup(group,subtreeL,familyT,strainNum2StrD,geneNames):
    '''Verbose print of a group.'''

    print("  Group",group.id)
    
    # get species nodes subtended by this mrca
    speciesNodesL=trees.leafList(subtreeL[group.mrca])

    # put everything in lists.
    printL=[]
    printL.append(['Family'])
    for node in speciesNodesL:
        printL[0].append(strainNum2StrD[node])
    for fam in group.familyL:
        newRow=[]
        newRow.append(str(fam))
        for node in speciesNodesL:
            ct,genesL = familyT[fam].famGeneT[node]
            newRow.append(",".join([geneNames.numToName(gene) for gene in genesL]))
        printL.append(newRow)
    printTable(printL,4)


def vPrintGroups(groupL,subtreeL,familyT,strainNum2StrD,geneNames):
    '''Print a list of groups.'''
    print("Summary of groups")
    printGroupLSummary(groupL)
    print("Print outs of each group")
    for group in groupL:
        vPrintGroup(group,subtreeL,familyT,strainNum2StrD,geneNames)
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

def createFam2GroupD(groupL):
    '''Given groupL, our list of groups, create a dictionary that maps
family number to group number.
    '''
    fam2GroupD={}
    for groupsAtNodeL in groupL:
        for group in groupsAtNodeL:
            for famNum in group.familyL:
                fam2GroupD[famNum]=group
    return fam2GroupD

def getGenesInFamily(family,subtreeL,familyT):
    '''Given a family, return a list of all genes in it (in numerical form).'''
    mrca = familyT[family].mrca
    leavesL=trees.leafList(subtreeL[mrca])

    genesL=[]
    for leaf in leavesL:

        geneT=familyT[family].famGeneT[leaf][1]
        genesL.extend(geneT)

    return genesL

    
def printScoreMatrix(family,subtreeL,familyT,geneNames,G):
    '''Print a matrix of scores between all the genes in a family. Scores
are provided by the graph G.'''

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
                row.append(format(data['score'],".3f"))
        rowsL.append(row)

    printTable(rowsL,2)

def printOutsideFamilyScores(family,subtreeL,familyT,geneNames,simG,synScoresG):
    '''Given a family, print scores for all non-family members with a
connection to genes in family. Scores are provided in the network
simG.
    '''
    
    familyGenesL = getGenesInFamily(family,subtreeL,familyT)

    rowL = []
    for gene in familyGenesL:
        geneName = geneNames.numToName(gene)
        for edge in simG.edges_iter(gene):
            if not edge[1] in familyGenesL:
                # this connection is with a gene outside the family
                otherGeneName = geneNames.numToName(edge[1])
                simData=simG.get_edge_data(gene,edge[1])
                synScoresData=synScoresG.get_edge_data(gene,edge[1])
                rowL.append([geneName,otherGeneName,format(simData['score'],".3f"),format(synScoresData['score'],".3f")])

    rowL.sort(key=lambda x: x[2],reverse=True) # sort by score
    rowL.insert(0,['----------','-----------','----------','-------'])
    rowL.insert(0,['Inside fam','Outside fam','Similarity','Synteny'])

                
    print("Printing all scores with non-family members")
    printTable(rowL,2)

                
def printFamNeighb(family,synWSize,subtreeL,familyT,geneOrderT,gene2FamD,fam2GroupD,geneDescriptionsD,geneNames,strainNum2StrD):
    '''Family Neighborhood. Given a family and an mrca for it, print out
the families nearby in each of its species and going either direction,
within synWSize of the last gene. If a dict of gene descriptions is
given, then include these in printout.
    '''

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
                    end = ind + synWSize
                    st = ind-synWSize if ind-synWSize>0 else 0 # st can't be less than 0
                    neighbGenesL=contig[st:end]
                    break
                except ValueError:
                    continue

            rowsL=[]
            for tempGene in neighbGenesL:
                geneName=geneNames.numToName(tempGene)
                famNum=gene2FamD[tempGene]
                group=fam2GroupD[famNum]
                descrip = geneDescriptionsD[geneName] if geneName in geneDescriptionsD else ''

                # if its the one in the family we're querying, mark with *
                if tempGene == gene:
                    geneName = '*'+geneName
                else:
                    geneName = ' '+geneName
                infoL = [geneName,"fam:"+str(famNum),"gr:"+str(group.id),"mrca:"+strainNum2StrD[group.mrca],descrip]

                rowsL.append(infoL)

            printTable(rowsL,4)


def printFam(synWSize,family):
    '''This is a wrapper to provide an easy way to print relevant info on
a family. For ease of use, we take only two arguments, assuming all
the other required stuff is available at the top level. Family is the
numerical identifier of a family. synWSize is the size of the genomic
window we will include when printing out genomic context in each
species.
    '''

    print()
    print("Matrix of similarity scores between genes in the family")
    printScoreMatrix(family,subtreeL,familyT,geneNames,simG)
    print()
    print()
        
    print("Matrix of synteny scores between genes in the family")
    printScoreMatrix(family,subtreeL,familyT,geneNames,synScoresG)
    print()
    print()
    
    printOutsideFamilyScores(family,subtreeL,familyT,geneNames,simG,synScoresG)
    print()
    print()

    print("Synteny information")
    printFamNeighb(family,synWSize,subtreeL,familyT,geneOrderT,gene2FamD,fam2GroupD,geneDescriptionsD,geneNames,strainNum2StrD)

def printGroupsAtNode(nodeStr):
    '''This is a wrapper to provide an easy way to print all the groups at
a particular node in the tree. For ease of use, we take only a node
number as argument, assuming all the other required stuff is available
at the top level.
    '''
    node = strainStr2NumD[nodeStr]
    vPrintGroups(groupByNodeL[node],subtreeL,familyT,strainNum2StrD,geneNames)

    
if __name__ == "__main__":

    paramFN=sys.argv[1]

    params = __import__(paramFN.replace('.py', ''))

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)
    
    # load groups
    groupByNodeL=readGroupOut(params.groupOutFN,tree,strainStr2NumD)
    
    # get familyT etc.
    geneNames = genomes.geneNames(params.geneOrderFN,strainStr2NumD,strainNum2StrD)

    
    geneDescriptionsD = genomes.createGeneDescriptionsD(params.geneDescriptionsFN)

    familyT = families.readFamilies(params.familyFN,tree,geneNames,strainStr2NumD)

    
    gene2FamD=createGene2FamD(familyT)
    fam2GroupD=createFam2GroupD(groupByNodeL)

    # subtree list
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()


    geneOrderT=genomes.createGeneOrderTs(params.geneOrderFN,geneNames,subtreeL,strainStr2NumD)

    # scores
    simG = scores.readGraph(params.scoresFN,geneNames)
    synScoresG = scores.readGraph(params.synScoresFN,geneNames)
    
