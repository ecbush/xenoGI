import sys
from multiprocessing import Pool
from . import trees
from . import genomes
from . import analysis
from .Island import *
from .Family import *
import math 

## Main function


def speedTestLocIsl(geneOrderT,geneNames,subtreeL,tree,paramD,familiesO,strainStr2NumD,strainNum2StrD):
    '''Test function to merge locus islands at a single node with a single thread.'''

    numThreads = paramD['numThreads']
    rootFocalClade = paramD['rootFocalClade']
    islandOutFN = paramD['islandOutFN']
    geneProximityRange = paramD['geneProximityRange']
    proximityThresholdMerge1 = paramD['proximityThresholdMerge1']
    rscThresholdMerge1 = paramD['rscThresholdMerge1']
    
    ## Temp parameters (later add to parameters.py)
    maxClusterSize = 50
    geneProximityRange = 2 # can adjust here for now. Should be >= 1.

    testNode = strainStr2NumD['i5']
    
    ##

    geneProximityD = genomes.createGeneProximityD(geneOrderT,geneProximityRange )
    locIslByNodeL=createLocIslByNodeL(familiesO,tree)
    numIslandsAtEachNodeAtStartL = [len(L) for L in locIslByNodeL]

    lFamNumL = []
    for liO in locIslByNodeL[testNode]:
        lFamNumL.extend(liO.locusFamilyL)

    print("num islands at node to start",len(locIslByNodeL[testNode]))
        
    subtree = subtreeL[testNode]
    
    _,cdD = costDiffDict((testNode,lFamNumL,familiesO,geneProximityD,proximityThresholdMerge1,subtree))

    # make famCostDiffL in a hacky way so only the one node is there
    famCostDiffL = [{} for i in range(trees.nodeCount(tree))]
    famCostDiffL[testNode] = cdD
    
    mergedL = mergeLocIslandsAtNode((locIslByNodeL[testNode],famCostDiffL,rscThresholdMerge1))

    print("num islands at node at end",len(mergedL))
    
    return mergedL
    
    
def makeLocusIslands(geneOrderT,geneNames,subtreeL,tree,paramD,familiesO,strainStr2NumD,strainNum2StrD, outputSummaryF):
    '''Parallelized wrapper to merge locus islands at different nodes.'''

    numThreads = paramD['numThreads']
    rootFocalClade = paramD['rootFocalClade']
    islandOutFN = paramD['islandOutFN']
    geneProximityRange = paramD['geneProximityRange']
    proximityThresholdMerge1 = paramD['proximityThresholdMerge1']
    rscThresholdMerge1 = paramD['rscThresholdMerge1']
    
    ## Temp parameters (later add to parameters.py)
    maxClusterSize = 50
    geneProximityRange = 2 # can adjust here for now. Should be >= 1.
    ##

    geneProximityD = genomes.createGeneProximityD(geneOrderT,geneProximityRange )
    locIslByNodeL=createLocIslByNodeL(familiesO,tree)
    numIslandsAtEachNodeAtStartL = [len(L) for L in locIslByNodeL]
    focalNodesL = getFocalNodesInOrderOfNumDescendants(tree,strainStr2NumD,rootFocalClade)

    famCostDiffMerge1L = createFamCostDiffL(tree,focalNodesL,locIslByNodeL,familiesO,geneProximityD,proximityThresholdMerge1,subtreeL,numThreads)

    
    print("num islands at start",numIslandsAtEachNodeAtStartL)


    ##  First merge (includes both merging in clusters and at mrca nodes)

    """

    # The clustering currently isn't working. So I'm going to skip it for now.

    locusIslandClusterL,singletonClusterL = createLocusIslandClusters(locIslByNodeL,focalNodesL,subtreeL,familiesO,geneProximityD,geneProximityRange,maxClusterSize)

    print("Num clusters:",len(locusIslandClusterL))
    print("Num singletons:",len(singletonClusterL))
    
    clusterLengths =  [len(L) for L in locusIslandClusterL]
    averageLen = sum(clusterLengths)/len(locusIslandClusterL)

    biggestC = max(clusterLengths)
    print("Biggest cluster has length of: ",biggestC)

    print("Average cluster size is: ", averageLen)

    # create argumentL to be passed to p.map and mergeLocIslandsAtNode
    argumentL = []
    for clusterL in locusIslandClusterL:
        argumentL.append((clusterL,geneProximityD,proximityThresholdMerge1,rscThresholdMerge1,subtreeL[clusterL[0].mrca],familiesO))
    p=Pool(numThreads)
    mergedL = p.map(mergeLocIslandsAtNode, argumentL) # run it

    # update locIslByNodeL with the merged nodes
    locIslByNodeL = updateIslandByNodeLEntries(locIslByNodeL,focalNodesL,mergedL)


    print("Cluster merge",[len(nodeL) for nodeL in locIslByNodeL])

    """

    ##  merge again at mrca nodes
    argumentL = []
    for mrcaNode in focalNodesL:
        argumentL.append((locIslByNodeL[mrcaNode],famCostDiffMerge1L,rscThresholdMerge1))
    p=Pool(numThreads)
    mergedL = p.map(mergeLocIslandsAtNode, argumentL) # run it

    # add the islands that were identified as singleton clusters
    #mergedL.append(singletonClusterL) # <-- Shouldn't this be extend?
    
    # again update locIslByNodeL
    locIslByNodeL = updateIslandByNodeLEntries(locIslByNodeL,focalNodesL,mergedL)

    numIslandsAtEachNodeAtEndL = [len(L) for L in locIslByNodeL]

    print("Node merge",numIslandsAtEachNodeAtEndL)
    
    ## Family improvement

    # not currently doing this
    
    ## Final merge

    # not currently doing this

    ## Summary and output
    
    # print summary of merging
    printSummary(focalNodesL,strainNum2StrD,numIslandsAtEachNodeAtStartL,numIslandsAtEachNodeAtEndL,outputSummaryF)
    
    # write islands
    writeIslands(locIslByNodeL,strainNum2StrD,islandOutFN)

    print("Islands written.",file=sys.stderr)

    return locIslByNodeL


def getFocalNodesInOrderOfNumDescendants(tree,strainStr2NumD,rootFocalClade):
    '''Get a list of all the nodes in the rootFocalClade. Then sort these
according to the number of leaves descending from them, biggest
first.'''

    focalSubtree = trees.subtree(tree,strainStr2NumD[rootFocalClade])
    focalNodesL = []
    for node in trees.nodeList(focalSubtree):
        # go through every node in the focal clade
        numDescend = trees.leafCount(trees.subtree(tree,node)) # count descendants
        focalNodesL.append((node,numDescend))
    focalNodesL.sort(key=lambda x: x[1],reverse=True) # sort on 2nd pos, biggest first
    # extract and return only the nodes
    return [node for node,numDescend in focalNodesL]

    
def updateIslandByNodeLEntries(locIslByNodeL,focalNodesL,mergedL):
    '''Given a list mergedL containing merged clusters of LocusIslands,
put back in locIslByNodeL. mergedL clusters all have mrca within
focalNodesL. We first blank the corresponding mrca location in
locIslByNodeL, then add the contents of mergedL.'''
    
    # blank out the entries in locIslByNodeL which are within the root focal clade
    for mrcaNode in focalNodesL:
        locIslByNodeL[mrcaNode] = []

    # fill them with LocusIslands from mergedL
    for clusterL in mergedL:
        locIslByNodeL[clusterL[0].mrca].extend(clusterL)

    return locIslByNodeL

## Cluster formation

def createLocusIslandClusters(locIslByNodeL,focalNodesL,subtreeL,familiesO,geneProximityD,proximityThreshold,maxClusterSize):
    '''For every node in the focal clade, take the set of single family
LocusIslands in locIslByNodeL. Break this up into smaller
clusters based on the chromosomal distances between members of the
familes. Each resulting cluster should contain LocusIslands that we
are likely to merge.'''

    locusIslandClusterL = [] #holds non-singleton clusters 
    singletonClusterL = [] #holds all of the islands from all MRCA nodes that are returned in clusters only containing one island 
    
    #loop through all of the mrcaNode lists 
    for islandsAtMrcaNodeLIndex in focalNodesL:
        
        islandsAtMrcaNodeL = locIslByNodeL[islandsAtMrcaNodeLIndex] 
        
        mrcaClusterList,mrcaSingletonClusterL = createMrcaNodeClusters(islandsAtMrcaNodeL,subtreeL, familiesO, geneProximityD,proximityThreshold, maxClusterSize,islandsAtMrcaNodeLIndex)

        singletonClusterL.extend(mrcaSingletonClusterL)
        
        locusIslandClusterL.extend(mrcaClusterList)
    
    print(type(mrcaClusterList[0]))
    print(type(mrcaSingletonClusterL[0]), "\n")
    return locusIslandClusterL,singletonClusterL

def createMrcaNodeClusters(islandsAtMrcaNodeL,subTreeL,familiesO,geneProximityD,proximityThreshold, maxClusterSize,islandsAtMrcaNodeLIndex):
    '''Takes in a Mrca node and other node information 
    and clusters the nodes within the mrca node that are likely 
    to merge'''

    mrcaClusterList = []
    mrcaSingletonClusterL = []
    
    while islandsAtMrcaNodeL:
        currentCluster = [] #create new empty cluster 
        firstIsland = islandsAtMrcaNodeL[0] #get the first island to go into the new cluster 
        
        currentCluster.append(firstIsland)
        islandsAtMrcaNodeL.remove(firstIsland)
        
        populateCluster(islandsAtMrcaNodeL,subTreeL,familiesO,geneProximityD,proximityThreshold, 
                    maxClusterSize,islandsAtMrcaNodeLIndex, currentCluster)
        
        if len(currentCluster) > 1:
            mrcaClusterList.append(currentCluster)
        else:
            mrcaSingletonClusterL.extend(currentCluster)

    return mrcaClusterList,mrcaSingletonClusterL

def populateCluster(islandsAtMrcaNodeL,subTreeL,familiesO,geneProximityD,proximityThreshold, 
                    maxClusterSize,islandsAtMrcaNodeLIndex, currentCluster):
    '''Takes in a island  and parses 
    though all of the nodes in the mrcaNode that is 
    not already used in another cluster if it is within 
    proximity adds it to newCluster'''
    
    index = 1
    while index<len(islandsAtMrcaNodeL):
        if len(currentCluster) > maxClusterSize:
            break 
            
        subTree = subTreeL[islandsAtMrcaNodeLIndex]

        currentFam = familiesO.getLocusFamily(currentCluster[0].locusFamilyL[0])

        comparisonFam = familiesO.getLocusFamily(islandsAtMrcaNodeL[index].locusFamilyL[0])

        proxB = proximitySubtree(currentFam,comparisonFam,geneProximityD,proximityThreshold,subTree)
           
        if proxB: 
            currentCluster.append(islandsAtMrcaNodeL[index])
            islandsAtMrcaNodeL.remove(islandsAtMrcaNodeL[index])
        else:
            index+=1
        
def proximitySubtree(lfam0,lfam1,geneProximityD,proximityThreshold,subtree):
    '''Given two gene families with the same mrcaNode, return boolean
indicating whether any of their genes are within proximityThreshold of
each other.
    '''
    if subtree[1] == ():
        return proximity(lfam0,lfam1,geneProximityD,proximityThreshold,subtree[0])
    else:
        left = proximitySubtree(lfam0,lfam1,geneProximityD,proximityThreshold,subtree[1])
        right = proximitySubtree(lfam0,lfam1,geneProximityD,proximityThreshold,subtree[2])
        return left or right
        

## Support functions

def createLocIslByNodeL(familiesO,tree):
    '''Create locus islands, one family each. Initially store islands
separately by mrca in a list of lists. (The index of the outer list
corresponds to the mrca)
    '''
    locIslByNodeL=[[] for i in range(trees.nodeCount(tree))]

    for lfO in familiesO.iterLocusFamilies():
        liO = LocusIsland(lfO.locusFamNum, lfO.lfMrca, [lfO.locusFamNum])
        locIslByNodeL[liO.mrca].append(liO)

    # sort each list by island number
    for i in range(len(locIslByNodeL)):
        locIslByNodeL[i].sort(key=lambda x: x.id)
    
    return locIslByNodeL


def createFamCostDiffL(tree,focalNodesL,locIslByNodeL,familiesO,geneProximityD,proximityThreshold,subtreeL,numThreads):
    '''Create a data structure which contrains the costDiff scores between
every pair of families at each mrca. This consists of a list, where
the index corresponds to mrca. The element at that index is a
dictionary of scores. It is keyed by a tuple of family numbers.'''

    argumentL = []
    for mrcaNode in focalNodesL:
        subtree = subtreeL[mrcaNode]
        lFamNumL = []
        for liO in locIslByNodeL[mrcaNode]:
            lFamNumL.extend(liO.locusFamilyL) # should only be one in the list at this point...
        
        argumentL.append((mrcaNode,lFamNumL,familiesO,geneProximityD,proximityThreshold,subtree))
        
    p=Pool(numThreads)
    mergedL = p.map(costDiffDict, argumentL) # run it

    famCostDiffL = [{} for i in range(trees.nodeCount(tree))]
    for mrcaNode,cdD in mergedL:
        famCostDiffL[mrcaNode] = cdD

    return famCostDiffL
    
def costDiffDict(argT):
    '''Create a dictionary of costDiff scores between all pairs of locus
families in lFamNumL.'''

    mrcaNode,lFamNumL,familiesO,geneProximityD,proximityThreshold,subtree=argT
    
    cdD={}
    for i in range(len(lFamNumL)-1):
        lf1 = familiesO.getLocusFamily(lFamNumL[i])
        for j in range(i+1,len(lFamNumL)):
            lf2 = familiesO.getLocusFamily(lFamNumL[j])
            cdsc = costDiff(lf1,lf2,geneProximityD,proximityThreshold,subtree)
            cdD[(lf1.locusFamNum,lf2.locusFamNum)] = cdsc
            cdD[(lf2.locusFamNum,lf1.locusFamNum)] = cdsc
    return mrcaNode,cdD
            
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
    
def rscore(li0,li1,famCostDiffL):
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
        cdD = famCostDiffL[li0.mrca] # costDiff score D at this mrca
        caseL = [-float('inf')]*4
        
        # case 0: last lfam in li0 vs. first lfam in li1
        caseL[0] = cdD[(li0.locusFamilyL[-1],li1.locusFamilyL[0])]

        # case 1: last lfam in li0 vs. last lfam in li1
        if len(li1.locusFamilyL) == 1:
            caseL[1] = caseL[0] # since first and last of li1 are same
        else:
            caseL[1] = cdD[(li0.locusFamilyL[-1],li1.locusFamilyL[-1])]
            
        # case 2: first lfam in li0 vs. first lfam in li1
        if len(li0.locusFamilyL) == 1:
            caseL[2] = caseL[0] # since first and last of li0 are same
        else:
            caseL[2] = cdD[(li0.locusFamilyL[0],li1.locusFamilyL[0])]
            
        # case 3: first lfam in li0 vs. last lfam in li1
        if len(li0.locusFamilyL) == 1 and len(li1.locusFamilyL) > 1:
            caseL[3] = caseL[1]
        elif len(li0.locusFamilyL) > 1 and len(li1.locusFamilyL) == 1:
            caseL[3] = caseL[2]
        elif len(li0.locusFamilyL) == 1 and len(li1.locusFamilyL) == 1:
            caseL[3] = caseL[0]
        else: # both longer than 1
            caseL[3] = cdD[(li0.locusFamilyL[0],li1.locusFamilyL[-1])]            
            
        return tuple(caseL)

def storeScore(li0,li1,famCostDiffL,scoreD):
    '''Calculate and store score in scoreD. We follow convention that key
in scoreD should always have the lower island id first.'''
    if li0.id < li1.id:
        key = li0.id,li1.id
        tempScoreT=rscore(li0,li1,famCostDiffL)
        # rscore returns different things depending on order
        # so we must be consistent with what we do in key.
    else:
        key = li1.id,li0.id
        tempScoreT=rscore(li1,li0,famCostDiffL)
        
    scoreD[key]=(max(tempScoreT),tempScoreT)

    
def createScoreD(locusIslandL,famCostDiffL):
    '''Create dictionary of scores between all locus islands at a single node
(and thus with the same mrca).'''
    scoreD={}
    for i in range(len(locusIslandL)-1):
        for j in range(i+1,len(locusIslandL)):
            storeScore(locusIslandL[i],locusIslandL[j],famCostDiffL,scoreD)
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

def addScores(scoreD,li0,locusIslandL,famCostDiffL):
    '''Get scores for locus island li0 against all other locus islands and add to
scoreD.'''
    for locIsl in locusIslandL:
        if locIsl.id != li0.id:
            storeScore(li0,locIsl,famCostDiffL,scoreD)

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
threshold. famCostDiffL contains precomputed costDiff scores between
families which we use to calculate rscore values between
islands. rscThreshold represents the threshold below which we no
longer merge islands.
    '''

    locusIslandL,famCostDiffL,rscThreshold = argT
    
    if len(locusIslandL) < 2:
        # nothing to merge
        return locusIslandL

    scoreD = createScoreD(locusIslandL,famCostDiffL)

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
        addScores(scoreD,li0,locusIslandL,famCostDiffL)

    return locusIslandL

## Output

def printSummary(focalNodesL,strainNum2StrD,numIslandsAtEachNodeAtStartL,numIslandsAtEachNodeAtEndL,outputSummaryF):
    '''Print a summary of the merging saying how many islands there were
at each node in the focal clade, before and after merging.'''

    print("Number of islands per node in focal clade: ",file=outputSummaryF)
    rowL=[]
    rowL.append(['Node','Before merge','After merge'])
    rowL.append(['----','------------','-----------'])
    for mrcaNode in focalNodesL:

        rowL.append([strainNum2StrD[mrcaNode],str(numIslandsAtEachNodeAtStartL[mrcaNode]),str(numIslandsAtEachNodeAtEndL[mrcaNode])])
    
    analysis.printTable(rowL,indent=2,fileF=outputSummaryF)

    return

def writeIslands(locIslByNodeL,strainNum2StrD,islandOutFN):
    '''Write the islands to a file'''
    f=open(islandOutFN,"w")
    for branch in range(len(locIslByNodeL)):
        for island in locIslByNodeL[branch]:
            print(island.fileStr(strainNum2StrD),file=f)
    f.close()


def readIslands(islandFN,tree,strainStr2NumD):
    '''Given a file name for a islands output file, load back
recreating locIslByNodeL.'''

    locIslByNodeL=[[] for i in range(trees.nodeCount(tree))]
    
    f=open(islandFN,'r')
    while True:
        s=f.readline()
        if s == '':
            break
        gr=str2Island(s.rstrip(),strainStr2NumD)
        locIslByNodeL[gr.mrca].append(gr)
    f.close()

    return locIslByNodeL
