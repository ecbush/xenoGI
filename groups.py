import sys
from multiprocessing import Pool
import trees,genomes
from Group import *
from Family import *


def createGroupL(familyT,tree):
    '''Greate groups, one family per group initially store groups
separately by mrca in a list of lists. (where the index of the outer
list corresponds to the mrca)'''
    groupL=[[] for i in range(trees.nodeCount(tree))]

    for fam in familyT:
        gr = Group(fam.id, fam.mrca, [fam.id])
        groupL[gr.mrca].append(gr)

    # sort each list by group number
    for i in range(len(groupL)):
        groupL[i].sort(key=lambda x: x.id)
    
    return groupL

            
## Distance calculations

def proximity(fam1,fam2,geneProximityD,proximityThreshold,node):
    '''Return True if any of the genes at node from fam1 are within
proximityThreshold of genes at that node for fam2. proximityThreshold
is measured in genes, e.g. 1 means the adjacent gene.
    '''
    for gn1 in fam1.famGeneT[node][1]:
        for gn2 in fam2.famGeneT[node][1]:
            if gn1<gn2: key = (gn1,gn2)
            else: key = (gn2,gn1)
            if key in geneProximityD and geneProximityD[key] <= proximityThreshold:
                return True
    return False
  
def rcost(fam1,fam2,geneProximityD,proximityThreshold,subtree,rootProximate):
    '''Parsimony based calculation of cost of evolutionary rearrangement
given fam1 and fam2 begin at root of subtree. Assume either proximate
(nearby) or not at root (specificed by rootProximate). Charge 1 for
each rearrangment.
    '''
    if subtree[1]==():
        prox=proximity(fam1,fam2,geneProximityD,proximityThreshold,subtree[0])
        if (prox and rootProximate) or (not prox and not rootProximate):
            # state given by rootProximate matches adjacency info in our data
            return 0
        else:
            return 1
    else:
        # left subtree
        left = rcost(fam1,fam2,geneProximityD,proximityThreshold,subtree[1],rootProximate)
        chLeft = 1 + rcost(fam1,fam2,geneProximityD,proximityThreshold,subtree[1],not rootProximate)

        # right subtree
        right = rcost(fam1,fam2,geneProximityD,proximityThreshold,subtree[2],rootProximate)
        chRight = 1 + rcost(fam1,fam2,geneProximityD,proximityThreshold,subtree[2],not rootProximate)

        return min(left,chLeft) + min(right,chRight)

def costDiff(fam1,fam2,geneProximityD,proximityThreshold,subtree):
    '''Given two families calculate the difference in rcost depending on
whether we assume the root is not proximate or proximate. fam1 and
fam2 are family tuples specifying the genes present in a family.
    '''
    t=rcost(fam1,fam2,geneProximityD,proximityThreshold,subtree,True)
    f=rcost(fam1,fam2,geneProximityD,proximityThreshold,subtree,False)
    return(f-t)
    
def rscore(g0,g1,geneProximityD,proximityThreshold,subtree,familyT):
    '''Returns the rearrangement score between groups g0 and g1. Considers
all four ways these groups could join (since there are two groups each
with two ends). This is the rcost given we start not proximate minus
the rcost given we start proximate.
    '''
    if g0.mrca != g1.mrca:
        return -float('inf')
    else:
        # we try combining in each of the 4 orientations (but if group
        # only has 1 family, then we can skip the reverse orientation
        # for that one.

        caseL = [-float('inf')]*4

        # case 0: last fam in g0 vs. first fam in g1
        caseL[0] = costDiff(familyT[g0.familyL[-1]],familyT[g1.familyL[0]],geneProximityD,proximityThreshold,subtree)

        # case 1: last fam in g0 vs. last fam in g1
        if len(g1.familyL) == 1:
            caseL[1] = caseL[0] # since first and last of g1 are same
        else:
            caseL[1] = costDiff(familyT[g0.familyL[-1]],familyT[g1.familyL[-1]],geneProximityD,proximityThreshold,subtree)

        # case 2: first fam in g0 vs. first fam in g1
        if len(g0.familyL) == 1:
            caseL[2] = caseL[0] # since first and last of g0 are same
        else:
            caseL[2] = costDiff(familyT[g0.familyL[0]],familyT[g1.familyL[0]],geneProximityD,proximityThreshold,subtree)

        # case 3: first fam in g0 vs. last fam in g1
        if len(g0.familyL) == 1 and len(g1.familyL) > 1:
            caseL[3] = caseL[1]
        elif len(g0.familyL) > 1 and len(g1.familyL) == 1:
            caseL[3] = caseL[2]
        elif len(g0.familyL) == 1 and len(g1.familyL) == 1:
            caseL[3] = caseL[0]
        else: # both longer than 1
            caseL[3] = costDiff(familyT[g0.familyL[0]],familyT[g1.familyL[-1]],geneProximityD,proximityThreshold,subtree)

        return tuple(caseL)

def storeScore(g0,g1,geneProximityD,proximityThreshold,subtree,scoreD,familyT):
    '''Calculate and store score in scoreD. We follow convention that key
in scoreD should always have the lower group id first.'''
    if g0.id < g1.id:
        key = g0.id,g1.id
        tempScoreT=rscore(g0,g1,geneProximityD,proximityThreshold,subtree,familyT)
        # rscore returns different things depending on order
        # so we must be consistent with what we do in key.
    else:
        key = g1.id,g0.id
        tempScoreT=rscore(g1,g0,geneProximityD,proximityThreshold,subtree,familyT)
        
    scoreD[key]=(max(tempScoreT),tempScoreT)

    
def createScoreD(groupL,geneProximityD,proximityThreshold,subtree,familyT):
    '''Create dictionary of scores between all groups at a single node
(and thus with the same mrca).'''
    scoreD={}
    for i in range(len(groupL)-1):
        for j in range(i+1,len(groupL)):
            storeScore(groupL[i],groupL[j],geneProximityD,proximityThreshold,subtree,scoreD,familyT)
    return scoreD


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

def addScores(scoreD,g0,groupsAtThisNodeL,geneProximityD,proximityThreshold,subtree,familyT):
    '''Get scores for group g0 against all other groups and add to scoreD.'''
    for gr in groupsAtThisNodeL:
        if gr.id != g0.id:
            storeScore(g0,gr,geneProximityD,proximityThreshold,subtree,scoreD,familyT)

def searchGroupsByID(listOfGroups,id):
    '''Search for a group with id equal to id. Return the index of the
first we find, and None if their isn't one.
    '''
    for i,gr in enumerate(listOfGroups):
        if gr.id == id: return i,gr
    return None,None

    
def mergeGroupsAtNode(argT):
    '''Given a list of groups at one node (groupL) iteratively merge until
there are no more pairwise scores above threshold. We consider pairs
of (proximity threshold, rscore threshold) which are specified in the
list proxThreshL inside argT. So we might first do all merging for
proximity threshold 1 (adjacent genes), merging until there are no
groups whose pairwise rscore is >= to rscore threshold. Then we
proceed to the next (proximity threshold, rscore threshold) pair,
which might have for example proximity threshold 2 (we consider genes
nearby that are separated by 1 gene).
    '''

    groupL,geneProximityD,proxThreshL,subtree,familyT = argT
    
    if len(groupL) < 2:
        # nothing to merge
        return groupL

    # we consider every pair of proximity threshold and rscore
    # threshold in proxThreshL
    for proximityThreshold,rscThreshold in proxThreshL:

        scoreD = createScoreD(groupL,geneProximityD,proximityThreshold,subtree,familyT)
        
        while True:
            sc,groupPairT,scoreT=maxScore(scoreD)
            if sc < rscThreshold:
                break

            ind0,g0 = searchGroupsByID(groupL,groupPairT[0])
            ind1,g1 = searchGroupsByID(groupL,groupPairT[1])

            g0.merge(g1,scoreT.index(sc))

            # delete g1
            del groupL[ind1]

            # remove all scores in scoreD that involve g0 or g1
            delScores(scoreD,g0.id,g1.id)

            # calculate new scores for g0 against all other groups
            addScores(scoreD,g0,groupL,geneProximityD,proximityThreshold,subtree,familyT)

    return groupL

## Output

def writeGroups(groupByNodeL,strainNum2StrD,groupOutFN):
    '''Write the groups to a file'''
    f=open(groupOutFN,"w")
    for branch in range(len(groupByNodeL)):
        for group in groupByNodeL[branch]:
            print(group.fileStr(strainNum2StrD),file=f)
    f.close()


def readGroups(groupFN,tree,strainStr2NumD):
    '''Given a file name for a groups output file, load back
recreating groupByNodeL.'''

    groupByNodeL=[[] for i in range(trees.nodeCount(tree))]
    
    f=open(groupFN,'r')
    while True:
        s=f.readline()
        if s == '':
            break
        gr=str2Group(s.rstrip(),strainStr2NumD)
        groupByNodeL[gr.mrca].append(gr)
    f.close()

    return groupByNodeL

    
## Main function

def makeGroups(geneOrderT,geneNames,subtreeL,tree,proxThreshL,familyT,numThreads,strainNum2StrD,groupOutFN, outputSummaryF):
    '''Parallelized wrapper to merge groups at different nodes.'''

    maxGeneProximityForGroup = max([thresh for thresh,rsc in proxThreshL])
    geneProximityD = genomes.createGeneProximityD(geneOrderT,maxGeneProximityForGroup)
    groupByNodeL=createGroupL(familyT,tree)

    print("Number of groups per node before merging: ", ' '.join([str(len(x)) for x in groupByNodeL]),file=outputSummaryF)


    ## create argumentL to be passed to p.map and mergeGroupsAtNode
    argumentL = []
    for mrcaNode in range(len(groupByNodeL)):
        argumentL.append((groupByNodeL[mrcaNode],geneProximityD,proxThreshL,subtreeL[mrcaNode],familyT))

    # run it
    p=Pool(numThreads)
    groupByNodeLMerged = p.map(mergeGroupsAtNode, argumentL) 

    #print("Did not merge core groups (last entries in groupByNodeL).",file=sys.stderr)
    print("Merging complete.",file=sys.stderr)
    print("Number of groups per node after merging: ", ' '.join([str(len(x)) for x in groupByNodeLMerged]),file=outputSummaryF)
    
    # write groups
    writeGroups(groupByNodeLMerged,strainNum2StrD,groupOutFN)

    print("Groups written.",file=sys.stderr)
