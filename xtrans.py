import sys
from tree import *
from genomes import *
from Group import *

   
## Load functions

# Families

def createFamilyStrainT(familyFN,tree,geneName2NumD,geneName2StrainNumD):
    '''Create a tuple representation of families. Input file consists of
two columns, first is family number, second is gene name (as a
string). We create a tuple, where the index is family number. The
value is another tuple (we'll refer to it as a 'family tuple'). This
family tuple has indexes corresponding to nodes on the tree. At each
position we have another tuple (gene count, [tuple of genes]).
    '''

    # load into raw lists
    rawFamL=[]
    geneL=[]
    f=open(familyFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.split()
        rawFamL.append(int(L[0]))
        geneL.append(L[1])
        
        
    # create initial version as list of lists
    tempFamiliesL=[]
    for i in range(max(rawFamL)+1):
        tempFamiliesL.append([[0,[]] for j in range(nodeCount(tree))])
        
    # now populate
    for i in range(len(rawFamL)):
        famNum=rawFamL[i]
        geneID=geneName2NumD[geneL[i]]
        strainNum=geneName2StrainNumD[geneL[i]]
        tempFamiliesL[famNum][strainNum][0]+=1
        tempFamiliesL[famNum][strainNum][1].append(geneID)
        
    # convert to tuple of tuples.
    newFamiliesL=[]
    for tempFam in tempFamiliesL:
        newTempFam=[]
        for famGeneCount,famGeneL in tempFam:
            newTempFam.append((famGeneCount,tuple(famGeneL)))
            
        newFamiliesL.append(tuple(newTempFam))

    return tuple(newFamiliesL)

## Group functions

def createGroupL(familyStrainT,tree):
    '''Greate groups, one family per group initially store groups separately by mrca in a list of lists.  (where the index of the outer list corresponds to the mrca)'''
    groupL=[[] for i in range(nodeCount(tree))]
    for i in range(1,len(familyStrainT)):
        # skip 0th because we have no 0 family number
        gr=Group(i,[i],familyStrainT=familyStrainT,tree=tree)
        groupL[gr.mrca].append(gr)

    # sort each list by group number
    for i in range(len(groupL)):
        groupL[i].sort(key=lambda x: x.id)
    
    return groupL

            
## Distance calculations

def isAdjacent(famT1,famT2,adjacencyS,node):
    '''Return True if any of the genes at node from fam1 are adjacent to
genes at that node for fam2.'''
    for gn1 in famT1[node][1]:
        for gn2 in famT2[node][1]:
            if gn1<gn2: key = (gn1,gn2)
            else: key = (gn2,gn1)
            if key in adjacencyS:
                return True
    return False
  
def rcost(famT1,famT2,adjacencyS,subtree,rootAdjacent):
    '''Parsimony based calculation of cost of evolutionary rearrangement
given fam1 and fam2 begin at root of subtree. Assume either adjacent
or not at root (specificed by rootAdjacent). Charge 1 for each
rearrangment.
    '''
    if subtree[1]==():
        isAdj=isAdjacent(famT1,famT2,adjacencyS,subtree[0])
        if (isAdj and rootAdjacent) or (not isAdj and not rootAdjacent):
            # state given by rootAdjacent matches adjacency info in our data
            return 0
        else:
            return 1
    else:
        # left subtree
        left = rcost(famT1,famT2,adjacencyS,subtree[1],rootAdjacent)
        chLeft = 1 + rcost(famT1,famT2,adjacencyS,subtree[1],not rootAdjacent)

        right = rcost(famT1,famT2,adjacencyS,subtree[2],rootAdjacent)
        chRight = 1 + rcost(famT1,famT2,adjacencyS,subtree[2],not rootAdjacent)

        return min(left,chLeft) + min(right,chRight)

def costDiff(famT1,famT2,adjacencyS,subtree):
    '''Given two families calculate the difference in rcost depending on
whether we assume the root is not adjacent or adjacent. famT1 and
famT2 are family tuples specifying the genes present in a family.
    '''
    t=rcost(famT1,famT2,adjacencyS,subtree,True)
    f=rcost(famT1,famT2,adjacencyS,subtree,False)
    return(f-t)


    
def score(g0,g1,adjacencyS,subtree):
    '''Returns the score between groups g0 and g1. Considers all four ways
these groups could join (since there are two groups each with two
ends). This is the rcost given we start not adjacent minus the rcost
given we start adjacent.
    '''
    if g0.mrca != g1.mrca:
        return -float('inf')
    else:
        # we try combining in each of the 4 orientations (but if group only has 1 family, then we can skip the reverse orientation for that one.

        caseL = [-float('inf')]*4

        # case 0: last fam in g0 vs. first fam in g1
        caseL[0] = costDiff(familyStrainT[g0.familyL[-1]],familyStrainT[g1.familyL[0]],adjacencyS,subtree)

        # case 1: last fam in g0 vs. last fam in g1
        if len(g1.familyL) == 1:
            caseL[1] = caseL[0] # since first and last of g1 are same
        else:
            caseL[1] = costDiff(familyStrainT[g0.familyL[-1]],familyStrainT[g1.familyL[-1]],adjacencyS,subtree)

        # case 2: first fam in g0 vs. first fam in g1
        if len(g0.familyL) == 1:
            caseL[2] = caseL[0] # since first and last of g0 are same
        else:
            caseL[2] = costDiff(familyStrainT[g0.familyL[0]],familyStrainT[g1.familyL[0]],adjacencyS,subtree)

        # case 3: first fam in g0 vs. last fam in g1
        if len(g0.familyL) == 1 and len(g1.familyL) > 1:
            caseL[3] = caseL[1]
        elif len(g0.familyL) > 1 and len(g1.familyL) == 1:
            caseL[3] = caseL[2]
        elif len(g0.familyL) == 1 and len(g1.familyL) == 1:
            caseL[3] = caseL[0]
        else: # both longer than 1
            caseL[3] = costDiff(familyStrainT[g0.familyL[0]],familyStrainT[g1.familyL[-1]],adjacencyS,subtree)

        return tuple(caseL)

def storeScore(g0,g1,adjacencyS,subtree,scoreD):
    '''Calculate and store score in scoreD. We follow convention that key
in scoreD should always have the lower group id first.'''
    if g0.id < g1.id:
        key = g0.id,g1.id
        tempScoreT=score(g0,g1,adjacencyS,subtree)
        # score returns different things depending on order
        # so we must be consisstent with what we do in key.
    else:
        key = g1.id,g0.id
        tempScoreT=score(g1,g0,adjacencyS,subtree)
    #if g1.id == 717:
    #    print(717,tempScoreT)
    #    print(g0)
    #    print(g1)
    scoreD[key]=(max(tempScoreT),tempScoreT)

    
def createScoreL(groupL,adjacencyS,subtreeL):
    '''Create list of scores, organized by mcra node. We only compare
groups with the same mrca. These scores are in a dictionary, which is
located at the index of the list corresponding to that mrca node.'''

    scoreL=[]
    for mrcaNode in range(len(groupL)-1): # -1 to skip core genes
        mrcaG=groupL[mrcaNode]
        scoreD={}
        for i in range(len(mrcaG)-1):
            for j in range(i+1,len(mrcaG)):
                if i == j: print("!")
                storeScore(mrcaG[i],mrcaG[j],adjacencyS,subtreeL[mrcaNode],scoreD)
                
        scoreL.append(scoreD)
    return scoreL


## Iterative merging

def maxScore(scoreD):
    '''Find the entry with the highest score in scoreD, and return the
corresponding key.'''
    bestSc=-float('inf')
    bestGroupPairT=()
    bestScoreT=()
    for groupPairT,(sc,scoreT) in scoreD.items():
        if sc > bestSc:
            bestSc=sc
            bestGroupPairT=groupPairT
            bestScoreT=scoreT
    return bestSc,bestGroupPairT,bestScoreT

def delScores(scoreD,g0ID,g1ID):
    '''Given two group ids from newly merged nodes, delete any entries in
score D that come from them.'''
    # find ones to delete
    toDelL=[]
    for key in scoreD.keys():
        if g0ID in key or g1ID in key:
            toDelL.append(key)

    # now delete them
    for key in toDelL:
        del scoreD[key]

def addScores(scoreD,g0,groupsAtThisNodeL,adjacencyS,subtree):
    '''Get scores for group g0 against all other groups and add to scoreD.'''
    for gr in groupsAtThisNodeL:
        if gr.id != g0.id:
            storeScore(g0,gr,adjacencyS,subtree,scoreD)

def searchGroupsByID(listOfGroups,id):
    '''Search for a group with id equal to id. Return the index of the
first we find, and None if their isn't one.
    '''
    for i,gr in enumerate(listOfGroups):
        if gr.id == id: return i,gr
    return None,None

    
def mergeGroupsAtNode(groupL,scoreL,adjacencyS,subtreeL,mrcaNode,threshold):
    '''Iteratively merge groups at mrcaNode.'''

    if len(groupL[mrcaNode]) < 2: return # nothing to merge
    
    sc=threshold

    iteration=0
    while True:
        #print("----iter",iteration)
        
        scoreD=scoreL[mrcaNode]
        sc,groupPairT,scoreT=maxScore(scoreD)

        if sc < threshold:
            break
        
        #print(sc,groupPairT,scoreT)

        ind0,g0 = searchGroupsByID(groupL[mrcaNode],groupPairT[0])
        ind1,g1 = searchGroupsByID(groupL[mrcaNode],groupPairT[1])

        #print("  group pair 0:",groupPairT[0],ind0,g0)
        #print("  group pair 1:",groupPairT[1],ind1,g1)
        
        g0.merge(g1,scoreT.index(sc))

        #print("merged group",g0)
        
        # delete g1
        del groupL[mrcaNode][ind1]

        
        # remove all scores in scoreL that involve g0 or g1
        delScores(scoreD,g0.id,g1.id)

        # calculate new scores for g0 against all other groups
        addScores(scoreD,g0,groupL[mrcaNode],adjacencyS,subtreeL[mrcaNode])

        iteration+=1

    return

## Output

def writeGroups(groupL,groupOutFN):
    '''Write the groups to a file'''
    f=open(groupOutFN,"w")
    for branch in range(len(groupL)):
        for group in groupL[branch]:
            print(group.fileStr(),file=f)



## Main

if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    # load data
    tree,strainStr2NumD,strainNum2StrD = readTree(params.treeFN)

    geneName2NumD,geneNum2NameD,geneName2StrainNumD = createGeneDs(params.geneOrderFN,strainStr2NumD)

    familyStrainT = createFamilyStrainT(params.familyFN,tree,geneName2NumD,geneName2StrainNumD)

    adjacencyS = createAdjacencySet(params.geneOrderFN,geneName2NumD)

    # create group list
    groupL=createGroupL(familyStrainT,tree)


    # subtree list
    subtreeL=createSubtreeL(tree)
    subtreeL.sort()


    ## cut off the last one, which will always be core families
    #groupL=groupL[:-1]

    print("Number of groups per node before merging: ", ' '.join([str(len(x)) for x in groupL]))
    
    # create score matrix
    print("Creating score matrix.",file=sys.stderr)
    scoreL=createScoreL(groupL,adjacencyS,subtreeL)

    
    # iteratively merge groups
    print("Begining merging.",file=sys.stderr)
    
    for i in range(len(groupL)-1):
        print("  Merging",i,file=sys.stderr)
        mergeGroupsAtNode(groupL,scoreL,adjacencyS,subtreeL,i,1)

    print("  Did not merge core groups (last entries in groupL).")
    print("Merging complete.",file=sys.stderr)

    print("Number of groups per node after merging: ", ' '.join([str(len(x)) for x in groupL]))

    # write groups
    writeGroups(groupL,params.groupOutFN)

    print("Groups written.",file=sys.stderr)

