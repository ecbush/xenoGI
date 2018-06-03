import sys
from multiprocessing import Pool
from . import trees
from . import genomes
from . import analysis
from .Island import *
from .Family import *

## Main function

def makeIslands(geneOrderT,geneNames,subtreeL,tree,proxThreshL,familyL,numThreads,strainStr2NumD,strainNum2StrD,rootFocalClade,islandOutFN, outputSummaryF):
    '''Parallelized wrapper to merge islands at different nodes.'''

    maxGeneProximityForIsland = max([thresh for thresh,rsc in proxThreshL])
    geneProximityD = genomes.createGeneProximityD(geneOrderT,maxGeneProximityForIsland)
    islandByNodeL=createIslandL(familyL,tree)

    focalSubtree = trees.subtree(tree,strainStr2NumD[rootFocalClade])
    focalNodesL=trees.nodeList(focalSubtree)
    

    ## create argumentL to be passed to p.map and mergeIslandsAtNode
    argumentL = []
    for mrcaNode in range(len(islandByNodeL)):
        if mrcaNode in focalNodesL:
            argumentL.append((islandByNodeL[mrcaNode],geneProximityD,proxThreshL,subtreeL[mrcaNode],familyL))
            
    # run it
    p=Pool(numThreads)
    islandByNodeLMerged = p.map(mergeIslandsAtNode, argumentL) 

    # islandByNodeLMerged has all islands from focal clade. Now add in
    # the single family islands from other nodes
    for mrcaNode in range(len(islandByNodeL)):
        if mrcaNode not in focalNodesL: 
            islandByNodeLMerged.append(islandByNodeL[mrcaNode])
            
    print("Merging complete.",file=sys.stderr)

    # print summary of merging
    printSummary(islandByNodeLMerged,focalNodesL,strainNum2StrD,islandByNodeL,outputSummaryF)
    
    # write islands
    writeIslands(islandByNodeLMerged,strainNum2StrD,islandOutFN)

    print("Islands written.",file=sys.stderr)

## Support functions

def createIslandL(familyL,tree):
    '''Greate islands, one family per island. Initially store islands
separately by mrca in a list of lists. (The index of the outer list
corresponds to the mrca)
    '''
    islandL=[[] for i in range(trees.nodeCount(tree))]

    for fam in familyL:
        isl = Island(fam.id, fam.mrca, [fam.id])
        islandL[isl.mrca].append(isl)

    # sort each list by island number
    for i in range(len(islandL)):
        islandL[i].sort(key=lambda x: x.id)
    
    return islandL

            
## Distance calculations

def proximity(fam1,fam2,geneProximityD,proximityThreshold,node):
    '''Return True if any of the genes at node from fam1 are within
proximityThreshold of genes at that node for fam2. proximityThreshold
is measured in genes, e.g. 1 means the adjacent gene.
    '''
    for gn1 in fam1.famGeneT[node]:
        for gn2 in fam2.famGeneT[node]:
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
    
def rscore(is0,is1,geneProximityD,proximityThreshold,subtree,familyL):
    '''Returns the rearrangement score between islands is0 and is1. Considers
all four ways these islands could join (since there are two islands each
with two ends). This is the rcost given we start not proximate minus
the rcost given we start proximate.
    '''
    if is0.mrca != is1.mrca:
        return -float('inf')
    else:
        # we try combining in each of the 4 orientations (but if island
        # only has 1 family, then we can skip the reverse orientation
        # for that one.

        caseL = [-float('inf')]*4

        # case 0: last fam in is0 vs. first fam in is1
        caseL[0] = costDiff(familyL[is0.familyL[-1]],familyL[is1.familyL[0]],geneProximityD,proximityThreshold,subtree)

        # case 1: last fam in is0 vs. last fam in is1
        if len(is1.familyL) == 1:
            caseL[1] = caseL[0] # since first and last of is1 are same
        else:
            caseL[1] = costDiff(familyL[is0.familyL[-1]],familyL[is1.familyL[-1]],geneProximityD,proximityThreshold,subtree)

        # case 2: first fam in is0 vs. first fam in is1
        if len(is0.familyL) == 1:
            caseL[2] = caseL[0] # since first and last of is0 are same
        else:
            caseL[2] = costDiff(familyL[is0.familyL[0]],familyL[is1.familyL[0]],geneProximityD,proximityThreshold,subtree)

        # case 3: first fam in is0 vs. last fam in is1
        if len(is0.familyL) == 1 and len(is1.familyL) > 1:
            caseL[3] = caseL[1]
        elif len(is0.familyL) > 1 and len(is1.familyL) == 1:
            caseL[3] = caseL[2]
        elif len(is0.familyL) == 1 and len(is1.familyL) == 1:
            caseL[3] = caseL[0]
        else: # both longer than 1
            caseL[3] = costDiff(familyL[is0.familyL[0]],familyL[is1.familyL[-1]],geneProximityD,proximityThreshold,subtree)

        return tuple(caseL)

def storeScore(is0,is1,geneProximityD,proximityThreshold,subtree,scoreD,familyL):
    '''Calculate and store score in scoreD. We follow convention that key
in scoreD should always have the lower island id first.'''
    if is0.id < is1.id:
        key = is0.id,is1.id
        tempScoreT=rscore(is0,is1,geneProximityD,proximityThreshold,subtree,familyL)
        # rscore returns different things depending on order
        # so we must be consistent with what we do in key.
    else:
        key = is1.id,is0.id
        tempScoreT=rscore(is1,is0,geneProximityD,proximityThreshold,subtree,familyL)
        
    scoreD[key]=(max(tempScoreT),tempScoreT)

    
def createScoreD(islandL,geneProximityD,proximityThreshold,subtree,familyL):
    '''Create dictionary of scores between all islands at a single node
(and thus with the same mrca).'''
    scoreD={}
    for i in range(len(islandL)-1):
        for j in range(i+1,len(islandL)):
            storeScore(islandL[i],islandL[j],geneProximityD,proximityThreshold,subtree,scoreD,familyL)
    return scoreD


## Iterative merging

def maxScore(scoreD):
    '''Find the entry with the highest score in scoreD, and return the
corresponding key.'''
    bestSc=-float('inf')
    bestIslandPairT=()
    bestScoreT=()
    for islandPairT,(sc,scoreT) in scoreD.items():
        if sc > bestSc:
            bestSc=sc
            bestIslandPairT=islandPairT
            bestScoreT=scoreT
    return bestSc,bestIslandPairT,bestScoreT

def delScores(scoreD,is0ID,is1ID):
    '''Given two island ids from newly merged nodes, delete any entries in
score D that come from them.'''
    # find ones to delete
    toDelL=[]
    for key in scoreD.keys():
        if is0ID in key or is1ID in key:
            toDelL.append(key)

    # now delete them
    for key in toDelL:
        del scoreD[key]

def addScores(scoreD,is0,islandsAtThisNodeL,geneProximityD,proximityThreshold,subtree,familyL):
    '''Get scores for island is0 against all other islands and add to
scoreD.'''
    for isl in islandsAtThisNodeL:
        if isl.id != is0.id:
            storeScore(is0,isl,geneProximityD,proximityThreshold,subtree,scoreD,familyL)

def searchIslandsByID(listOfIslands,id):
    '''Search for an island with id equal to id. Return the index of the
first we find, and None if their isn't one.
    '''
    for i,isl in enumerate(listOfIslands):
        if isl.id == id: return i,isl
    return None,None
    
def mergeIslandsAtNode(argT):
    '''Given a list of islands at one node (islandL) iteratively merge until
there are no more pairwise scores above threshold. We consider pairs
of (proximity threshold, rscore threshold) which are specified in the
list proxThreshL inside argT. So we might first do all merging for
proximity threshold 1 (adjacent genes), merging until there are no
islands whose pairwise rscore is >= to rscore threshold. Then we
proceed to the next (proximity threshold, rscore threshold) pair,
which might have for example proximity threshold 2 (we consider genes
nearby that are separated by 1 gene).
    '''

    islandL,geneProximityD,proxThreshL,subtree,familyL = argT
    
    if len(islandL) < 2:
        # nothing to merge
        return islandL

    # we consider every pair of proximity threshold and rscore
    # threshold in proxThreshL
    for proximityThreshold,rscThreshold in proxThreshL:

        scoreD = createScoreD(islandL,geneProximityD,proximityThreshold,subtree,familyL)
        
        while True:
            sc,islandPairT,scoreT=maxScore(scoreD)
            if sc < rscThreshold:
                break

            ind0,is0 = searchIslandsByID(islandL,islandPairT[0])
            ind1,is1 = searchIslandsByID(islandL,islandPairT[1])

            is0.merge(is1,scoreT.index(sc))

            # delete is1
            del islandL[ind1]

            # remove all scores in scoreD that involve is0 or is1
            delScores(scoreD,is0.id,is1.id)

            # calculate new scores for is0 against all other islands
            addScores(scoreD,is0,islandL,geneProximityD,proximityThreshold,subtree,familyL)

    return islandL

## Output

def printSummary(islandByNodeLMerged,focalNodesL,strainNum2StrD,islandByNodeL,outputSummaryF):
    '''Print a summary of the merging saying how many islands there were
at each node in the focal clade, before and after merging.'''

    focalPrintD = {}
    for islandGroup in islandByNodeLMerged:
        if islandGroup != []:
            mrca = islandGroup[0].mrca
            if mrca in focalNodesL:
                focalPrintD[mrca] = [strainNum2StrD[mrca],str(len(islandByNodeL[mrca])),str(len(islandGroup))]

    print("Number of islands per node in focal clade: ",file=outputSummaryF)
    rowL=[]
    rowL.append(['Node','Before merge','After merge'])
    rowL.append(['----','------------','-----------'])
    for node in focalNodesL:
        if node in focalPrintD:
            rowL.append(focalPrintD[node])
        else:
            rowL.append([strainNum2StrD[node],'0','0'])

    analysis.printTable(rowL,indent=2,fileF=outputSummaryF)

    return


def writeIslands(islandByNodeL,strainNum2StrD,islandOutFN):
    '''Write the islands to a file'''
    f=open(islandOutFN,"w")
    for branch in range(len(islandByNodeL)):
        for island in islandByNodeL[branch]:
            print(island.fileStr(strainNum2StrD),file=f)
    f.close()


def readIslands(islandFN,tree,strainStr2NumD):
    '''Given a file name for a islands output file, load back
recreating islandByNodeL.'''

    islandByNodeL=[[] for i in range(trees.nodeCount(tree))]
    
    f=open(islandFN,'r')
    while True:
        s=f.readline()
        if s == '':
            break
        gr=str2Island(s.rstrip(),strainStr2NumD)
        islandByNodeL[gr.mrca].append(gr)
    f.close()

    return islandByNodeL
