import sys
from Group import *
from xtrans import *

## Load groups

def readGroupOut(fn,tree):
    '''Given a file name for an htrans groups output file, load back
recreating groupL.'''

    groupL=[[] for i in range(nodeCount(tree))]
    
    f=open(fn,'r')
    while True:
        s=f.readline()
        if s == '':
            break
        gr=str2Group(s.rstrip())
        groupL[gr.mrca].append(gr)
    f.close()

    return groupL

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
        print(" "*indent + " | ".join(row))

        
def vPrintGroup(group,subtreeL,familyStrainT,strainNum2StrD,geneNum2NameD):
    '''Verbose print of a group.'''

    print("Group",group.id)
    
    # get species nodes subtended by this mrca
    speciesNodesL=leafList(subtreeL[group.mrca])

    # put everything in lists.
    printL=[]
    printL.append(['Family'])
    for node in speciesNodesL:
        printL[0].append(strainNum2StrD[node])
    for fam in group.familyL:
        newRow=[]
        newRow.append(str(fam))
        for node in speciesNodesL:
            ct,genesL = familyStrainT[fam][node]
            newRow.append(",".join([geneNum2NameD[gene] for gene in genesL]))
        printL.append(newRow)
    printTable(printL,2)


def vPrintGroups(groupL,subtreeL,familyStrainT,strainNum2StrD,geneNum2NameD):
    '''Print a list of groups.'''
    for group in groupL:
        vPrintGroup(group,subtreeL,familyStrainT,strainNum2StrD,geneNum2NameD)
        print('---')

def createGeneOrderTs(geneOrderFN,geneName2NumD,subtreeL,strainStr2NumD):
    '''Go though gene order file and get orderings into a set of tuples.'''
    f = open(geneOrderFN,'r')
    geneOrderL=[None for x in range(nodeCount(subtreeL[-1]))] # an index for each node
    while True:
        s = f.readline()
        if s == '':
            break
        s=s.rstrip()
        # note, our gene order format has contigs separated by \t, and
        # genes within them separated by a space character.
        L=s.split('\t')
        strain = L[0]
        gnNmContigL=[]
        for contig in L[1:]:
            geneNumT=tuple((geneName2NumD[g] for g in contig.split(' ')))
            gnNmContigL.append(geneNumT)
            
        geneOrderL[strainStr2NumD[strain]]=tuple(gnNmContigL)
    return tuple(geneOrderL)

def createGene2FamD(familyStrainT):
    '''Given the family information in familyStrainT, create a dictionary
gene2FamD which maps from gene number to family number.'''
    gene2FamD={}
    for famNum in range(len(familyStrainT)):
        for gnCt,geneT in familyStrainT[famNum]:
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
        
def printFamNeighb(family,mrca,wsize,subtreeL,familyStrainT,gene2FamD,fam2GroupD):
    '''Family Neighborhood. Given a family and an mrca for it, print out
the families nearby in each of its species and going either direction,
within wsize of the last gene.'''

    leavesL=leafList(subtreeL[mrca])

    for leaf in leavesL:

        print("Neighbors for family",family,"in",strainNum2StrD[leaf])

        geneT=familyStrainT[family][leaf][1]
        #if geneT[0] > 1:
        #    print("Right now can't handle cases where more than one gene in a species")
        #    return

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
                try:
                    geneName=geneNum2NameD[tempGene]
                    if tempGene == gene: # if its the one in the family we're querying, mark with *
                        geneName = '*'+geneName
                    else:
                        geneName = ' '+geneName
                    famNum=gene2FamD[tempGene]
                    group=fam2GroupD[famNum]
                    rowsL.append([geneName,str(famNum),str(group)])
                except KeyError:
                    # There are genes which are present in the
                    # ordering files, but not in the sequences
                    # files. This can lead to key errors here.
                    rowsL.append([' geneNum '+str(tempGene),'no family','no group'])
                    
            printTable(rowsL,2)


        
if __name__ == "__main__":

    paramFN=sys.argv[1]

    params = __import__(paramFN.replace('.py', ''))

    tree,strainStr2NumD,strainNum2StrD = readTree(params.treeFN)
    
    # load groups
    groupL=readGroupOut(params.groupOutFN,tree)

    
    # get familyStrainT etc.
    

    geneName2NumD,geneNum2NameD,geneName2StrainNumD = createGeneDs(params.geneOrderFN,strainStr2NumD)

    familyStrainT = createFamilyStrainT(params.familyFN,tree,geneName2NumD,geneName2StrainNumD)

    gene2FamD=createGene2FamD(familyStrainT)
    fam2GroupD=createFam2GroupD(groupL)

    # subtree list
    subtreeL=createSubtreeL(tree)
    subtreeL.sort()


    geneOrderT=createGeneOrderTs(params.geneOrderFN,geneName2NumD,subtreeL,strainStr2NumD)

    # Go.
    #vPrintGroups(groupL[2],subtreeL,familyStrainT,strainNum2StrD,geneNum2NameD)
    #printFamNeighb(3518,2,5,subtreeL,familyStrainT,gene2FamD,fam2GroupD)
    #printFamNeighb(610,2,5,subtreeL,familyStrainT,gene2FamD,fam2GroupD)
