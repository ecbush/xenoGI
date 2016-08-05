import sys
import genomes,trees,groups,scores
from Family import *
from Group import *

## Load groups

def readGroupOut(fn,tree):
    '''Given a file name for an htrans groups output file, load back
recreating groupByNodeL.'''

    groupByNodeL=[[] for i in range(trees.nodeCount(tree))]
    
    f=open(fn,'r')
    while True:
        s=f.readline()
        if s == '':
            break
        gr=str2Group(s.rstrip())
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
            
    
    
def vPrintGroup(group,subtreeL,familyStrainT,strainNum2StrD,geneNum2NameD):
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
            ct,genesL = familyStrainT[fam].famGeneT[node]
            newRow.append(",".join([geneNum2NameD[gene] for gene in genesL]))
        printL.append(newRow)
    printTable(printL,4)


def vPrintGroups(groupL,subtreeL,familyStrainT,strainNum2StrD,geneNum2NameD):
    '''Print a list of groups.'''
    print("Summary of groups")
    printGroupLSummary(groupL)
    print("Print outs of each group")
    for group in groupL:
        vPrintGroup(group,subtreeL,familyStrainT,strainNum2StrD,geneNum2NameD)
        print('  ---')

def createGene2FamD(familyStrainT):
    '''Given the family information in familyStrainT, create a dictionary
gene2FamD which maps from gene number to family number.'''
    gene2FamD={}
    for famNum in range(len(familyStrainT)):
        for gnCt,geneT in familyStrainT[famNum].famGeneT:
            for gene in geneT:
                gene2FamD[gene]=famNum
    return gene2FamD

def createFam2GroupD(groupL):
    '''Given groupL, our list of groups, create a dictionary that takes
maps family number to group number.'''
    fam2GroupD={}
    for groupsAtNodeL in groupL:
        for group in groupsAtNodeL:
            for famNum in group.familyL:
                fam2GroupD[famNum]=group
    return fam2GroupD

def printScoreMatrix(family,subtreeL,familyStrainT,geneNum2NameD,G):
    '''--p'''
    mrca = familyStrainT[family].mrca
    leavesL=trees.leafList(subtreeL[mrca])

    genesL=[]
    for leaf in leavesL:

        geneT=familyStrainT[family].famGeneT[leaf][1]
        genesL.extend(geneT)

    rowsL = []
    geneNamesL = [geneNum2NameD[gn] for gn in genesL]
    rowsL.append([''] + geneNamesL)
    
    for rowi,gn1 in enumerate(genesL):
        row = [geneNum2NameD[genesL[rowi]]]
        for gn2 in genesL:
            data=G.get_edge_data(gn1,gn2)
            if data == None:
                row.append('-')
            else:
                row.append(format(data['score'],".3f"))
        rowsL.append(row)

    print("Printing family",family)
    printTable(rowsL,2)



    
def printFamNeighb(family,wsize,subtreeL,familyStrainT,geneOrderT,gene2FamD,fam2GroupD,geneDescriptionsD):
    '''Family Neighborhood. Given a family and an mrca for it, print out
the families nearby in each of its species and going either direction,
within wsize of the last gene. If a dict of gene descriptions is
given, then include these in printout.
    '''

    mrca = familyStrainT[family].mrca
    leavesL=trees.leafList(subtreeL[mrca])

    for leaf in leavesL:

        print("Neighbors for family",family,"in",strainNum2StrD[leaf])

        geneT=familyStrainT[family].famGeneT[leaf][1]
        if geneT == ():
            print("  There are no family members in this species.")

        for gene in geneT:

            print("  For",geneNum2NameD[gene])
            
            neighbGenesL=[]
            for contig in geneOrderT[leaf]:
                try:
                    ind=contig.index(gene)
                    neighbGenesL=contig[ind-wsize:ind+wsize]
                    break
                except ValueError:
                    continue

            rowsL=[]
            for tempGene in neighbGenesL:

                if tempGene not in geneNum2NameD:
                    # There are genes which are present in the
                    # ordering files, but not in the sequence
                    # files. This can lead to key errors here.
                    infoL = [' geneNum '+str(tempGene),'no family','no group','no mrca','no descrip']
                else:
                    geneName=geneNum2NameD[tempGene]
                    famNum=gene2FamD[tempGene]
                    group=fam2GroupD[famNum]
                    descrip = geneDescriptionsD[geneName] if geneName in geneDescriptionsD else ''

                    # if its the one in the family we're querying, mark with *
                    if tempGene == gene:
                        geneName = '*'+geneName
                    else:
                        geneName = ' '+geneName
                    infoL = [geneName,"fam:"+str(famNum),"gr:"+str(group.id),"mrca:"+str(group.mrca),descrip]

                rowsL.append(infoL)

                    
            printTable(rowsL,2)


        
if __name__ == "__main__":

    paramFN=sys.argv[1]

    params = __import__(paramFN.replace('.py', ''))

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)
    
    # load groups
    groupByNodeL=readGroupOut(params.groupOutFN,tree)

    
    # get familyStrainT etc.
    

    geneName2NumD,geneNum2NameD,geneName2StrainNumD = genomes.createGeneDs(params.geneOrderFN,strainStr2NumD)
    geneDescriptionsD = genomes.createGeneDescriptionsD(params.geneDescriptionsFN)

    familyStrainT = groups.createFamilyStrainT(params.familyFN,tree,geneName2NumD,geneName2StrainNumD)

    gene2FamD=createGene2FamD(familyStrainT)
    fam2GroupD=createFam2GroupD(groupByNodeL)

    # subtree list
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()


    geneOrderT=genomes.createGeneOrderTs(params.geneOrderFN,geneName2NumD,subtreeL,strainStr2NumD)

    # scores
    simG = scores.createSimilarityGraph(params.scoresFN,geneName2NumD)
    synScoresG = scores.createSimilarityGraph(params.synScoresFN,geneName2NumD)

    
    # Go.
    #vPrintGroups(groupByNodeL[2],subtreeL,familyStrainT,strainNum2StrD,geneNum2NameD)
    #printFamNeighb(2569,5,subtreeL,familyStrainT,geneOrderT,gene2FamD,fam2GroupD,geneDescriptionsD = geneDescriptionsD)

