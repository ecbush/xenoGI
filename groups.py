import trees

from Group import *
from Family import *

# Families

def createFamilyStrainT(familyFN,tree,geneName2NumD,geneName2StrainNumD):
    '''Create a tuple representation of families. Input file consists of
one family per line. Family number, mcra node number, and genes in
family. We create a tuple, where the index is family number. The value
is an object of class Family.
    '''

    rawFamL=[]
    maxFamNum=0
    f=open(familyFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.split()
        famNum=int(L[0])
        mrca = int(L[1])
        genesL = L[2:]

        rawFamL.append(Family(famNum,mrca,genesL,trees.nodeCount(tree),geneName2NumD,geneName2StrainNumD))

        if famNum > maxFamNum:
            maxFamNum = famNum

    # ensure each family is put at index corresponding to its fam num
    famL = [None for i in range(maxFamNum+1)]
    for fam in rawFamL:
        famL[fam.id] = fam
    
    return tuple(famL)


## Group functions

def createGroupL(familyStrainT,tree):
    '''Greate groups, one family per group initially store groups
separately by mrca in a list of lists.  (where the index of the outer
list corresponds to the mrca)'''
    groupL=[[] for i in range(trees.nodeCount(tree))]

    for fam in familyStrainT:
        gr = Group(fam.id, fam.mrca, [fam.id])
        groupL[gr.mrca].append(gr)

    # sort each list by group number
    for i in range(len(groupL)):
        groupL[i].sort(key=lambda x: x.id)
    
    return groupL

            
## Distance calculations

def isAdjacent(fam1,fam2,adjacencyS,node):
    '''Return True if any of the genes at node from fam1 are adjacent to
genes at that node for fam2.'''
    for gn1 in fam1.famGeneT[node][1]:
        for gn2 in fam2.famGeneT[node][1]:
            if gn1<gn2: key = (gn1,gn2)
            else: key = (gn2,gn1)
            if key in adjacencyS:
                return True
    return False
  
def rcost(fam1,fam2,adjacencyS,subtree,rootAdjacent):
    '''Parsimony based calculation of cost of evolutionary rearrangement
given fam1 and fam2 begin at root of subtree. Assume either adjacent
or not at root (specificed by rootAdjacent). Charge 1 for each
rearrangment.
    '''
    if subtree[1]==():
        isAdj=isAdjacent(fam1,fam2,adjacencyS,subtree[0])
        if (isAdj and rootAdjacent) or (not isAdj and not rootAdjacent):
            # state given by rootAdjacent matches adjacency info in our data
            return 0
        else:
            return 1
    else:
        # left subtree
        left = rcost(fam1,fam2,adjacencyS,subtree[1],rootAdjacent)
        chLeft = 1 + rcost(fam1,fam2,adjacencyS,subtree[1],not rootAdjacent)

        right = rcost(fam1,fam2,adjacencyS,subtree[2],rootAdjacent)
        chRight = 1 + rcost(fam1,fam2,adjacencyS,subtree[2],not rootAdjacent)

        return min(left,chLeft) + min(right,chRight)

def costDiff(fam1,fam2,adjacencyS,subtree):
    '''Given two families calculate the difference in rcost depending on
whether we assume the root is not adjacent or adjacent. fam1 and
fam2 are family tuples specifying the genes present in a family.
    '''
    t=rcost(fam1,fam2,adjacencyS,subtree,True)
    f=rcost(fam1,fam2,adjacencyS,subtree,False)
    return(f-t)


    
def score(g0,g1,adjacencyS,subtree,familyStrainT):
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

def storeScore(g0,g1,adjacencyS,subtree,scoreD,familyStrainT):
    '''Calculate and store score in scoreD. We follow convention that key
in scoreD should always have the lower group id first.'''
    if g0.id < g1.id:
        key = g0.id,g1.id
        tempScoreT=score(g0,g1,adjacencyS,subtree,familyStrainT)
        # score returns different things depending on order
        # so we must be consisstent with what we do in key.
    else:
        key = g1.id,g0.id
        tempScoreT=score(g1,g0,adjacencyS,subtree,familyStrainT)
    #if g1.id == 717:
    #    print(717,tempScoreT)
    #    print(g0)
    #    print(g1)
    scoreD[key]=(max(tempScoreT),tempScoreT)

    
def createScoreL(groupL,adjacencyS,subtreeL,familyStrainT):
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
                storeScore(mrcaG[i],mrcaG[j],adjacencyS,subtreeL[mrcaNode],scoreD,familyStrainT)
                
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

def addScores(scoreD,g0,groupsAtThisNodeL,adjacencyS,subtree,familyStrainT):
    '''Get scores for group g0 against all other groups and add to scoreD.'''
    for gr in groupsAtThisNodeL:
        if gr.id != g0.id:
            storeScore(g0,gr,adjacencyS,subtree,scoreD,familyStrainT)

def searchGroupsByID(listOfGroups,id):
    '''Search for a group with id equal to id. Return the index of the
first we find, and None if their isn't one.
    '''
    for i,gr in enumerate(listOfGroups):
        if gr.id == id: return i,gr
    return None,None

    
def mergeGroupsAtNode(groupL,scoreL,adjacencyS,subtreeL,mrcaNode,threshold,familyStrainT):
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
        addScores(scoreD,g0,groupL[mrcaNode],adjacencyS,subtreeL[mrcaNode],familyStrainT)

        iteration+=1

    return

## Output

def writeGroups(groupL,groupOutFN):
    '''Write the groups to a file'''
    f=open(groupOutFN,"w")
    for branch in range(len(groupL)):
        for group in groupL[branch]:
            print(group.fileStr(),file=f)



