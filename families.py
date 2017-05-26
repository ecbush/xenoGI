# Functions for a modiefied version of the PhiGs algorithm
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# (We've added the use of synteny)
import sys
import trees
from Family import *


## Main function

def families(tree,subtreeL,geneNames,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,synAdjustThresh,synAdjustExtent,familyFN,strainNum2StrD, outputSummaryF):
    '''Given a graph of genes and their similarity scores find families
using a PhiGs-like algorithm, with synteny also considered.
    '''

    nodeOrderL=createNodeProcessOrderList(tree)

    # initialize scoresO.nodeConnectL for use below
    scoresO.createNodeConnectL(geneNames)
    
    geneUsedL = [False for x in geneNames.nums]
    multiGeneFamilyL=[]
    
    for node,lnode,rnode in nodeOrderL:
        if lnode != None:
            # not a tip
            subtree=subtreeL[node]
            leftS,rightS,outgroupS = createLRSets(subtree,geneNames)
            seedL = createSeedL(leftS,rightS,scoresO,minNormThresh,minCoreSynThresh,minSynThresh)
            for seed in seedL:
                seedSimScore,g1,g2 = seed

                if seedSimScore == -float('inf'):
                    # we've gotten to the point in the seed list with
                    # genes having no match on the other branch
                    break
                else:
                    family,possibleErrorCt=getFamily(seedSimScore,g1,g2,scoresO,leftS,rightS,minNormThresh,minCoreSynThresh,minSynThresh,synAdjustThresh,synAdjustExtent)
                    
                    if any((geneUsedL[gene] for gene in family)):
                        continue
                    else:
                        # none of the genes in this family used yet
                        for gene in family:
                            geneUsedL[gene] = True
                        multiGeneFamilyL.append((node,family,possibleErrorCt))

    # Create Family objects 
    genesInMultiGeneFamsS = set()
    familyL=[]
    famNum = 0
    for mrca,famS,possibleErrorCt in multiGeneFamilyL:
        genesInMultiGeneFamsS.update(famS) # add all genes in fam
        familyL.append(Family(famNum,mrca,famS,trees.nodeCount(tree),geneNames,possibleErrorCt ))
        famNum+=1

    multiGeneFamNum=famNum
    print("Number of multigene families",multiGeneFamNum,file=outputSummaryF)

    # Add on remaining genes as single gene families.
    for gene in geneNames.nums: 
        if not gene in genesInMultiGeneFamsS:
            mrca = geneNames.numToStrainNum(gene)
            familyL.append(Family(famNum,mrca,[gene],trees.nodeCount(tree),geneNames,None))
            famNum+=1

    print("Number of single gene families",famNum-multiGeneFamNum,file=outputSummaryF)
    print("Number of total families",famNum,file=outputSummaryF)
    
    writeFamilies(familyL,geneNames,strainNum2StrD,familyFN)

    return tuple(familyL)

## Support functions

def createNodeProcessOrderList(tree):
    '''Given a tree, output a list specifying nodes in pre-order (ancestors
before their descendants). For each node we give a tuple (node #, left
node #, right node #).
    '''
    if tree[1] == ():
        return [(tree[0], None, None)]
    else:
        l = createNodeProcessOrderList(tree[1])
        r = createNodeProcessOrderList(tree[2])
        return [(tree[0], tree[1][0], tree[2][0])] + l + r


def createLRSets(tree,geneNames):
    '''For every gene in our data, put it into one of three sets. Left,
right, or outgroup. Genes in the left set are found in a species on
the left branch of tree.'''

    leftSpeciesS=set(trees.leafList(tree[1]))
    rightSpeciesS=set(trees.leafList(tree[2]))
    
    leftS=set()
    rightS=set()
    outgroupS=set()
    for geneNum in geneNames.nums: # all genes
        strain=geneNames.numToStrainNum(geneNum)
        if strain in leftSpeciesS:
            leftS.add(geneNum)
        elif strain in rightSpeciesS:
            rightS.add(geneNum)
        else:
            outgroupS.add(geneNum)

    return(leftS,rightS,outgroupS)

def closestMatch(gene,S,scoresO,minNormThresh,minCoreSynThresh,minSynThresh):
    '''Find the closest match to gene among the genes in the set S in the
graph scoresO. Eliminate any matches that have a norm score below
minNormThresh, a coreSynSc below minCoreSynThresh, or a synteny score
below minSynThresh.
    '''
    bestGene=None
    bestEdgeScore = -float('inf')
    for otherGene in scoresO.getConnectionsGene(gene):
        if otherGene in S:
            if scoresO.getScoreByEndNodes(gene,otherGene,'rawSc') > bestEdgeScore:
                if scoresO.getScoreByEndNodes(gene,otherGene,'normSc') > minNormThresh:
                    if scoresO.getScoreByEndNodes(gene,otherGene,'coreSynSc') > minCoreSynThresh:
                        if scoresO.getScoreByEndNodes(gene,otherGene,'synSc') > minSynThresh:
                            bestEdgeScore = scoresO.getScoreByEndNodes(gene,otherGene,'rawSc')
                            bestGene = otherGene
    return bestEdgeScore, gene, bestGene
    
def createSeedL(leftS,rightS,scoresO,minNormThresh,minCoreSynThresh,minSynThresh):
    '''Create a list which has the closest match for each gene on the
opposite side of the tree. e.g. if a gene is in tree[1] then we're
looking for the gene in tree[2] with the closest match. We eliminate
any matches that are below threshold for normalized or syntenty
scores. For each gene we get this closest match, put in a list, sort,
and return.
    '''
    seedL=[]
    for gene in leftS:
        seedL.append(closestMatch(gene,rightS,scoresO,minNormThresh,minCoreSynThresh,minSynThresh))
    for gene in rightS:
        seedL.append(closestMatch(gene,leftS,scoresO,minNormThresh,minCoreSynThresh,minSynThresh))
    seedL.sort(reverse=True)
    return seedL

def getFamily(seedSimScore,g1,g2,scoresO,leftS,rightS,minNormThresh,minCoreSynThresh,minSynThresh,synAdjustThresh,synAdjustExtent):
    '''Based on a seed (seedScore, g1, g2) search for a family. Using the
PhiGs approach, we collect all genes which are closer to members of
the family than the two seeds are from each other. We have a
normScore, synScore, and coreSynScore thresholds below which we will
not add genes. We also have a synteny adjustment of the raw score
where we make the raw score between a pair a bit better if their
synteny is above synAdjustThresh. In general, if a gene has syntenic
connectsions to genes already in the family, this makes us more
confident that this gene belongs in the family. Returns a set
containing genes in the family.
    '''
    alreadySearchedS = set()
    notYetSearchedS = set([g1,g2])

    possibleErrorCt = 0
    while len(notYetSearchedS) > 0:
        newGenesS=set()
        for famGene in notYetSearchedS:
            for newGene in scoresO.getConnectionsGene(famGene):
                if newGene in leftS or newGene in rightS:
                    # it is from a species descended from the node
                    # we're working on
                    if newGene not in alreadySearchedS and newGene not in newGenesS:
                        addIt,err = addGene(famGene,newGene,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,seedSimScore,synAdjustThresh,synAdjustExtent)
                        possibleErrorCt += err
                        if addIt:
                            newGenesS.add(newGene)
                    
            alreadySearchedS.add(famGene)
        notYetSearchedS = newGenesS
    return alreadySearchedS,possibleErrorCt


def addGene(famGene,newGene,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,seedSimScore,synAdjustThresh,synAdjustExtent):
    '''Given famGene that is inside a family, and a newGene we are
considering adding, check the various scores to determine if we should
add it. Return boolean.
    '''

    # get scores
    rawSc = scoresO.getScoreByEndNodes(famGene,newGene,'rawSc')
    normSc = scoresO.getScoreByEndNodes(famGene,newGene,'normSc')
    synSc = scoresO.getScoreByEndNodes(famGene,newGene,'synSc')
    coreSynSc = scoresO.getScoreByEndNodes(famGene,newGene,'coreSynSc')


    # determine whether or not to add newGene
    # norm theshold is a sort of sanity check. Is the level of
    # similarity just too low?
    normB = normSc >= minNormThresh
    
    # two types of synteny threshold, core synteny and regular synteny
    coreB = coreSynSc >= minCoreSynThresh
    synB = synSc >= minSynThresh
    
    # compare raw score with seed score, adjusting for synteny
    if synSc > synAdjustThresh:
        # its above the syn score adjustment
        # threshold, so increase rawSc a bit. This
        # addresses a problem with closely related
        # families where the seed score is very
        # similar. Sometimes by chance things
        # that should have been added weren't
        # because they weren't more similar than
        # an already very similar seed.
        # (Modification of PhiGs)
        rawAdj = rawSc * synAdjustExtent
        if rawAdj > 1: rawAdj = 1 # truncate back to 1
    else:
        rawAdj = rawSc
    rawB = rawAdj >= seedSimScore
        
    addIt = normB and coreB and synB and rawB

    # determine if it's a possible error and return magic numbers
    # here...but we'll be replacing this function with a more
    # systematic way run after family formation.
    if isPossibleError(addIt,normSc,minNormThresh,1):
        return addIt,1
    elif isPossibleError(addIt,coreSynSc,minCoreSynThresh,.05):
        return addIt,1
    elif isPossibleError(addIt,synSc,minSynThresh,1):
        return addIt,1
    elif isPossibleError(addIt,rawAdj,seedSimScore,.1):
        return addIt,1
    else: return addIt,0
    
def isPossibleError(addIt,sc,threshold,errorScoreDetermIncrement):
    '''Determine if this was a near miss. If we added it, did we almost
    not do so? If we didn't add it, did we almost add it? If it's a
    near miss, it will contribute to our error score.
    '''
    if addIt:
        if not (sc - errorScoreDetermIncrement) >= threshold:
            # It was over, but now is below threshold. Thus its a
            # possible error
            return True
    else:
        if (sc + errorScoreDetermIncrement) >= threshold:
            # It was below, but now its over threshold. Thus its's a
            # possible error.
            return True
    return False

    
def addGeneOnceAndFuture(famGene,newGene,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,seedSimScore,synAdjustThresh,synAdjustExtent):
    '''Given famGene that is inside a family, and newGene we are
considering adding, check the various scores to determine if we should
add it. Return boolean.
    '''
    addIt = False
    
    rawSc = scoresO.getScoreByEndNodes(famGene,newGene,'rawSc')
    normSc = scoresO.getScoreByEndNodes(famGene,newGene,'normSc')
    synSc = scoresO.getScoreByEndNodes(famGene,newGene,'synSc')
    coreSynSc = scoresO.getScoreByEndNodes(famGene,newGene,'coreSynSc')

    if normSc < minNormThresh:
        # doesn't meet score threshold for family
        # formation. (Modification of PhiGs)
        pass
    elif coreSynSc < minCoreSynThresh:
        # doesn't meet min core synteny requirement for
        # family formation. (Modification of PhiGs)
        pass
    elif synSc < minSynThresh:
        # doesn't meet min synteny requirement for
        # family formation. (Modification of PhiGs)
        pass
    elif rawSc >= seedSimScore:
        # If its within the seed distance, add it
        # (basic PhiGs approach). we have the =
        # there in case seedSimScore is 1.
        addIt = True
    elif synSc > synAdjustThresh:
        # its above the syn score adjustment
        # threshold, so increase rawSc a bit. This
        # addresses a problem with closely related
        # families where the seed score is very
        # similar. Sometimes by chance things
        # that should have been added weren't
        # because they weren't more similar than
        # an already very similar seed.
        # (Modification of PhiGs)
        adjSc = rawSc * synAdjustExtent
        if adjSc > 1: adjSc = 1 # truncate back to 1
        if adjSc >= seedSimScore:
            addIt = True
    return addIt




## Input/output

def writeFamilies(familyL,geneNames,strainNum2StrD,fileName):
    '''Write all gene families to fileName, one family per line.'''
    f=open(fileName,'w')
    for fam in familyL:
        f.write(fam.fileStr(geneNames,strainNum2StrD)+'\n')
    f.close()


def readFamilies(familyFN,tree,geneNames,strainStr2NumD):
    '''Read the family file named familyFN, creating a tuple of family
objects where the index corresponds to family number.'''
    
    famL=[]
    maxFamNum=0
    f=open(familyFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.split()
        famNum=int(L[0])
        mrca = strainStr2NumD[L[1]]
        possibleErrorCt = L[2]
        if possibleErrorCt == 'None':
            possibleErrorCt = None
        else:
            possibleErrorCt = int(possibleErrorCt)
        
        genesL = L[3:]

        famL.append(Family(famNum,mrca,genesL,trees.nodeCount(tree),geneNames,possibleErrorCt))
    f.close()
    return tuple(famL)
