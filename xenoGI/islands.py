import sys
from multiprocessing import Pool
from . import trees
from . import genomes
from . import analysis
from .Island import *
from .Family import *

## Main function

def makeLocusIslands(geneOrderT,geneNames,subtreeL,tree,paramD,familiesO,strainStr2NumD,strainNum2StrD, outputSummaryF):
    '''Parallelized wrapper to merge locus islands at different nodes.'''

    proxThreshL = paramD['proxThreshL']
    numThreads = paramD['numThreads']
    rootFocalClade = paramD['rootFocalClade']
    islandOutFN = paramD['islandOutFN']
    
    maxGeneProximityForIsland = max([thresh for thresh,rsc in proxThreshL])
    geneProximityD = genomes.createGeneProximityD(geneOrderT,maxGeneProximityForIsland)
    locusIslandByNodeL=createLocusIslandL(familiesO,tree)

    focalSubtree = trees.subtree(tree,strainStr2NumD[rootFocalClade])
    focalNodesL=trees.nodeList(focalSubtree)
    

    ## create argumentL to be passed to p.map and mergeLocIslandsAtNode
    argumentL = []
    for mrcaNode in range(len(locusIslandByNodeL)):
        if mrcaNode in focalNodesL:
            argumentL.append((locusIslandByNodeL[mrcaNode],geneProximityD,proxThreshL,subtreeL[mrcaNode],familiesO))
            
    # run it
    p=Pool(numThreads)
    locusIslandByNodeLMerged = p.map(mergeLocIslandsAtNode, argumentL) 

    # locusIslandByNodeLMerged has all islands from focal clade. Now add in
    # the single family islands from other nodes
    for mrcaNode in range(len(locusIslandByNodeL)):
        if mrcaNode not in focalNodesL: 
            locusIslandByNodeLMerged.append(locusIslandByNodeL[mrcaNode])

    # print summary of merging
    printSummary(locusIslandByNodeLMerged,focalNodesL,strainNum2StrD,locusIslandByNodeL,outputSummaryF)

    # write islands
    writeIslands(locusIslandByNodeLMerged,strainNum2StrD,islandOutFN)

    print("Islands written.",file=sys.stderr)

    return locusIslandByNodeLMerged

## Support functions

def createLocusIslandL(familiesO,tree):
    '''Create locus islands, one family each. Initially store islands
separately by mrca in a list of lists. (The index of the outer list
corresponds to the mrca)
    '''
    locusIslandL=[[] for i in range(trees.nodeCount(tree))]

    for lfO in familiesO.iterLocusFamilies():
        liO = LocusIsland(lfO.locusFamNum, lfO.lfMrca, [lfO.locusFamNum])
        locusIslandL[liO.mrca].append(liO)

    # sort each list by island number
    for i in range(len(locusIslandL)):
        locusIslandL[i].sort(key=lambda x: x.id)
    
    return locusIslandL

            
## Distance calculations

def proximity(lfam1,lfam2,geneProximityD,proximityThreshold,strain):
    '''Return True if any of the genes at strain from lfam1 are within
proximityThreshold of genes at that strain for lfam2. proximityThreshold
is measured in genes, e.g. 1 means the adjacent gene.
    '''
    for gn1 in lfam1.iterGenesByStrain(strain):
        for gn2 in lfam2.iterGenesByStrain(strain):
            if gn1<gn2: key = (gn1,gn2)
            else: key = (gn2,gn1)
            if key in geneProximityD and geneProximityD[key] <= proximityThreshold:
                return True
    return False
  
def rcost(lfam1,lfam2,geneProximityD,proximityThreshold,subtree,rootProximate):
    '''Parsimony based calculation of cost of evolutionary rearrangement
given lfam1 and lfam2 begin at root of subtree. Assume either proximate
(nearby) or not at root (specificed by rootProximate). Charge 1 for
each rearrangment.
    '''
    if subtree[1]==():
        prox=proximity(lfam1,lfam2,geneProximityD,proximityThreshold,subtree[0])
        if (prox and rootProximate) or (not prox and not rootProximate):
            # state given by rootProximate matches adjacency info in our data
            return 0
        else:
            return 1
    else:
        # left subtree
        left = rcost(lfam1,lfam2,geneProximityD,proximityThreshold,subtree[1],rootProximate)
        chLeft = 1 + rcost(lfam1,lfam2,geneProximityD,proximityThreshold,subtree[1],not rootProximate)

        # right subtree
        right = rcost(lfam1,lfam2,geneProximityD,proximityThreshold,subtree[2],rootProximate)
        chRight = 1 + rcost(lfam1,lfam2,geneProximityD,proximityThreshold,subtree[2],not rootProximate)

        return min(left,chLeft) + min(right,chRight)

def costDiff(lfam1,lfam2,geneProximityD,proximityThreshold,subtree):
    '''Given two families calculate the difference in rcost depending on
whether we assume the root is not proximate or proximate. lfam1 and
lfam2 are family tuples specifying the genes present in a family.
    '''
    t=rcost(lfam1,lfam2,geneProximityD,proximityThreshold,subtree,True)
    f=rcost(lfam1,lfam2,geneProximityD,proximityThreshold,subtree,False)
    return(f-t)
    
def rscore(li0,li1,geneProximityD,proximityThreshold,subtree,familiesO):
    '''Returns the rearrangement score between LocusIslands li0 and
li1. Considers all four ways these locus islands could join (since
there are two locus islands each with two ends). This is the rcost
given we start not proximate minus the rcost given we start proximate.
    '''
    if li0.mrca != li1.mrca:
        return -float('inf')
    else:
        # we try combining in each of the 4 orientations (but if locus
        # island only has 1 family, then we can skip the reverse
        # orientation for that one.

        caseL = [-float('inf')]*4
        
        # case 0: last lfam in li0 vs. first lfam in li1
        caseL[0] = costDiff(familiesO.getLocusFamily(li0.locusFamilyL[-1]),familiesO.getLocusFamily(li1.locusFamilyL[0]),geneProximityD,proximityThreshold,subtree)

        # case 1: last lfam in li0 vs. last lfam in li1
        if len(li1.locusFamilyL) == 1:
            caseL[1] = caseL[0] # since first and last of li1 are same
        else:
            caseL[1] = costDiff(familiesO.getLocusFamily(li0.locusFamilyL[-1]),familiesO.getLocusFamily(li1.locusFamilyL[-1]),geneProximityD,proximityThreshold,subtree)

        # case 2: first lfam in li0 vs. first lfam in li1
        if len(li0.locusFamilyL) == 1:
            caseL[2] = caseL[0] # since first and last of li0 are same
        else:
            caseL[2] = costDiff(familiesO.getLocusFamily(li0.locusFamilyL[0]),familiesO.getLocusFamily(li1.locusFamilyL[0]),geneProximityD,proximityThreshold,subtree)

        # case 3: first lfam in li0 vs. last lfam in li1
        if len(li0.locusFamilyL) == 1 and len(li1.locusFamilyL) > 1:
            caseL[3] = caseL[1]
        elif len(li0.locusFamilyL) > 1 and len(li1.locusFamilyL) == 1:
            caseL[3] = caseL[2]
        elif len(li0.locusFamilyL) == 1 and len(li1.locusFamilyL) == 1:
            caseL[3] = caseL[0]
        else: # both longer than 1
            caseL[3] = costDiff(familiesO.getLocusFamily(li0.locusFamilyL[0]),familiesO.getLocusFamily(li1.locusFamilyL[-1]),geneProximityD,proximityThreshold,subtree)
                                
        return tuple(caseL)

def storeScore(li0,li1,geneProximityD,proximityThreshold,subtree,scoreD,familiesO):
    '''Calculate and store score in scoreD. We follow convention that key
in scoreD should always have the lower island id first.'''
    if li0.id < li1.id:
        key = li0.id,li1.id
        tempScoreT=rscore(li0,li1,geneProximityD,proximityThreshold,subtree,familiesO)
        # rscore returns different things depending on order
        # so we must be consistent with what we do in key.
    else:
        key = li1.id,li0.id
        tempScoreT=rscore(li1,li0,geneProximityD,proximityThreshold,subtree,familiesO)
        
    scoreD[key]=(max(tempScoreT),tempScoreT)

    
def createScoreD(locusIslandL,geneProximityD,proximityThreshold,subtree,familiesO):
    '''Create dictionary of scores between all locus islands at a single node
(and thus with the same mrca).'''
    scoreD={}
    for i in range(len(locusIslandL)-1):
        for j in range(i+1,len(locusIslandL)):
            storeScore(locusIslandL[i],locusIslandL[j],geneProximityD,proximityThreshold,subtree,scoreD,familiesO)
    return scoreD


## Iterative merging

def maxScore(scoreD):
    '''Find the entry with the highest score in scoreD, and return the
corresponding key.'''
    bestSc=-float('inf')
    bestLocIslandPairT=()
    bestScoreT=()
    for locIslandPairT,(sc,scoreT) in scoreD.items():
        if sc > bestSc:
            bestSc=sc
            bestLocIslandPairT=locIslandPairT
            bestScoreT=scoreT
    return bestSc,bestLocIslandPairT,bestScoreT

def delScores(scoreD,li0ID,li1ID):
    '''Given two locus island ids from newly merged nodes, delete any entries in
score D that come from them.'''
    # find ones to delete
    toDelL=[]
    for key in scoreD.keys():
        if li0ID in key or li1ID in key:
            toDelL.append(key)

    # now delete them
    for key in toDelL:
        del scoreD[key]

def addScores(scoreD,li0,islandsAtThisNodeL,geneProximityD,proximityThreshold,subtree,familyL):
    '''Get scores for locus island li0 against all other locus islands and add to
scoreD.'''
    for locIsl in islandsAtThisNodeL:
        if locIsl.id != li0.id:
            storeScore(li0,locIsl,geneProximityD,proximityThreshold,subtree,scoreD,familyL)

def searchLocIslandsByID(listOfLocIslands,id):
    '''Search for a locus island with id equal to id. Return the index of the
first we find, and None if their isn't one.
    '''
    for i,isl in enumerate(listOfLocIslands):
        if isl.id == id: return i,isl
    return None,None
    
def mergeLocIslandsAtNode(argT):
    '''Given a list of locus islands at one node (locusIslandL)
iteratively merge until there are no more pairwise scores above
threshold. We consider pairs of (proximity threshold, rscore
threshold) which are specified in the list proxThreshL inside argT. So
we might first do all merging for proximity threshold 1 (adjacent
genes), merging until there are no islands whose pairwise rscore is >=
to rscore threshold. Then we proceed to the next (proximity threshold,
rscore threshold) pair, which might have for example proximity
threshold 2 (we consider genes nearby that are separated by 1 gene).
    '''

    locusIslandL,geneProximityD,proxThreshL,subtree,familiesO = argT
    
    if len(locusIslandL) < 2:
        # nothing to merge
        return locusIslandL

    # we consider every pair of proximity threshold and rscore
    # threshold in proxThreshL
    for proximityThreshold,rscThreshold in proxThreshL:

        scoreD = createScoreD(locusIslandL,geneProximityD,proximityThreshold,subtree,familiesO)
        
        while True:
            sc,islandPairT,scoreT=maxScore(scoreD)
            if sc < rscThreshold:
                break

            ind0,li0 = searchLocIslandsByID(locusIslandL,islandPairT[0])
            ind1,li1 = searchLocIslandsByID(locusIslandL,islandPairT[1])

            li0.merge(li1,scoreT.index(sc))

            # delete li1
            del locusIslandL[ind1]

            # remove all scores in scoreD that involve li0 or li1
            delScores(scoreD,li0.id,li1.id)

            # calculate new scores for li0 against all other islands
            addScores(scoreD,li0,locusIslandL,geneProximityD,proximityThreshold,subtree,familiesO)

    return locusIslandL

## Output

def printSummary(locusIslandByNodeLMerged,focalNodesL,strainNum2StrD,locusIslandByNodeL,outputSummaryF):
    '''Print a summary of the merging saying how many islands there were
at each node in the focal clade, before and after merging.'''

    focalPrintD = {}
    for islandGroup in locusIslandByNodeLMerged:
        if islandGroup != []:
            mrca = islandGroup[0].mrca
            if mrca in focalNodesL:
                focalPrintD[mrca] = [strainNum2StrD[mrca],str(len(locusIslandByNodeL[mrca])),str(len(islandGroup))]

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


def writeIslands(locusIslandByNodeL,strainNum2StrD,islandOutFN):
    '''Write the islands to a file'''
    f=open(islandOutFN,"w")
    for branch in range(len(locusIslandByNodeL)):
        for island in locusIslandByNodeL[branch]:
            print(island.fileStr(strainNum2StrD),file=f)
    f.close()


def readIslands(islandFN,tree,strainStr2NumD):
    '''Given a file name for a islands output file, load back
recreating locusIslandByNodeL.'''

    locusIslandByNodeL=[[] for i in range(trees.nodeCount(tree))]
    
    f=open(islandFN,'r')
    while True:
        s=f.readline()
        if s == '':
            break
        gr=str2Island(s.rstrip(),strainStr2NumD)
        locusIslandByNodeL[gr.mrca].append(gr)
    f.close()

    return locusIslandByNodeL
