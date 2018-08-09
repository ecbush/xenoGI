# Functions for a modified version of the PhiGs algorithm
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# (We've added the use of synteny)
import sys
from . import trees,scores
from .Family import *
from .analysis import printTable

## Main function

def createFamiliesO(tree,strainNum2StrD,scoresO,geneNames,aabrhL,singleStrainFamilyThresholdAdjust,subtreeL,minNormThresh,minCoreSynThresh,minSynThresh,synAdjustThresh,synAdjustExtent,outputSummaryF,familyFN):
    '''Given a graph of genes and their similarity scores find families
using a PhiGs-like algorithm, with synteny also considered.

Here's a hierarchy of (major) function calls below this:
createFamiliesO
  createAllFamiliesDescendingFromInternalNode
    createFamilyFromSeed
  createAllFamiliesAtTip
  createAllLocusFamiliesDescendingFromInternalNode
     createAllLocusFamiliesAtOneInternalNode
        createLocusFamilyFromSeed
     createAllLocusFamiliesAtTip
    '''

    # initialize scoresO.nodeConnectL and scoresO.ScoreSummaryD for
    # use below
    strainNumsL=sorted(leaf for leaf in trees.leafList(tree))
    scoresO.createNodeConnectL(geneNames)
    scoresO.createAabrhScoreSummaryD(strainNumsL,aabrhL,geneNames)

    # create an object of class Families to store this in.
    familiesO = Families(tree)
    famNumCounter = 0
    locusFamNumCounter = 0

    # other assorted things we'll need
    geneUsedL = [False for x in geneNames.nums]
    nodeGenesL = createNodeGenesL(tree,geneNames) # has genes divided by node
    tipFamilyRawThresholdD = getTipFamilyRawThresholdD(tree,scoresO,singleStrainFamilyThresholdAdjust)
    
    for familyMrca,lchild,rchild in createNodeProcessOrderList(tree):
        # this is preorder, so we get internal nodes before tips
        if lchild != None:
            # not a tip
                        
            geneUsedL,locusFamNumCounter,famNumCounter,familiesO = createAllFamiliesDescendingFromInternalNode(subtreeL,familyMrca,nodeGenesL,scoresO,geneNames,minNormThresh,synAdjustThresh,synAdjustExtent,geneUsedL,familiesO,famNumCounter,locusFamNumCounter,minCoreSynThresh,minSynThresh)

        else:
            geneUsedL,locusFamNumCounter,famNumCounter,familiesO = createAllFamiliesAtTip(nodeGenesL,familyMrca,geneUsedL,tipFamilyRawThresholdD,scoresO,geneNames,minNormThresh,synAdjustThresh,synAdjustExtent,familiesO,minCoreSynThresh,minSynThresh,famNumCounter,locusFamNumCounter)

       
    # Write family formation summary file

    summaryL=[]
    summaryL.append(["Total number of families",str(len(familiesO.familiesD))])
    summaryL.append(["Total number of LocusFamilies",str(len(familiesO.locusFamiliesD))])

    singleLfFams = 0
    multipleLfFams = 0
    singleGeneFams = 0
    multipleGeneFams = 0
    for fam in familiesO.iterFamilies():
        # locus families
        if len(fam.getLocusFamilies()) == 1:
            singleLfFams += 1
        else:
            multipleLfFams += 1
        # genes
        if len(fam.getGeneNums()) == 1:
            singleGeneFams += 1
        else:
            multipleGeneFams += 1


    summaryL.append(["Number of families with one LocusFamily",str(singleLfFams)])
    summaryL.append(["Number of families with multiple LocusFamilies",str(multipleLfFams)])
    summaryL.append(["Number of families with one gene",str(singleGeneFams)])
    summaryL.append(["Number of families with multiple genes",str(multipleGeneFams)])
    
    printTable(summaryL,indent=0,fileF=outputSummaryF)
    
    writeFamilies(familiesO,geneNames,strainNum2StrD,familyFN)

    return familiesO

## Support functions

def createNodeGenesL(tree,geneNames):
    '''Create a data structure to organize genes by strain. Returns a list
where the indices correspond to strain number and the elements are
sets.
    '''
    nodeGenesL=[set() for i in range(trees.nodeCount(tree))]
    for gene in geneNames.nums:
        strain = geneNames.numToStrainNum(gene)
        nodeGenesL[strain].add(gene)
    return nodeGenesL

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

def getTipFamilyRawThresholdD(tree,scoresO,singleStrainFamilyThresholdAdjust):
    '''Return a dictionary containing a raw score threshold for each
tip. This threshold is for use in forming families on the tip,
defining the minimum distance within which we will combine two gene
into a family.'''

    tipFamilyRawThresholdD = {}
    for leaf in trees.leafList(tree):
        # get average score at core genes for neighbors

        # put in call to get average...
        threshold,std = scores.getNearestNeighborAverageScore(leaf,tree,scoresO)
        
        # multiply in an adjustment parameter (since average core gene
        # scores of neighbors would be way too high)
        threshold *= singleStrainFamilyThresholdAdjust

        tipFamilyRawThresholdD[leaf] = threshold

    return tipFamilyRawThresholdD


def createAllFamiliesDescendingFromInternalNode(subtreeL,familyMrca,nodeGenesL,scoresO,geneNames,minNormThresh,synAdjustThresh,synAdjustExtent,geneUsedL,familiesO,famNumCounter,locusFamNumCounter,minCoreSynThresh,minSynThresh):
    '''Creates all Families and subsidiary LocusFamilies descending from
the node rooted familyMrca. Basic parts of the Phigs algorithm are
here. Creating the seeds, and using them to get a family. (With very
minor changes in the use of the norm score and synteny
adjustment). But then we divide the family into LocusFamilies, which
is not Phigs.
    '''

    subtree=subtreeL[familyMrca]
    
    # for use in createAllLocusFamiliesDescendingFromInternalNode call below
    familySubtreeNodeOrderL = createNodeProcessOrderList(subtree)

    leftS,rightS = createLRSets(subtreeL,familyMrca,nodeGenesL,None)

    seedL = createSeedL(leftS,rightS,scoresO,geneNames,minNormThresh)
    for seed in seedL:
        # each seed corresponds to a prospective gene family.
        seedRawSc,seedG1,seedG2 = seed

        if seedRawSc == -float('inf'):
            # we've gotten to the point in the seed list with
            # genes having no match on the other branch
            break
        else:
            # getting initial family, using only raw score and synteny bump
            famS=createFamilyFromSeed(seedG1,seedG2,geneUsedL,scoresO,leftS,rightS,geneNames,seedRawSc,minNormThresh,synAdjustThresh,synAdjustExtent)

            if famS == None:
                # one of the genes the the family was already
                # used, so createFamilyFromSeed returned None
                continue
            else:
                # none of the genes in famS used yet
                for gene in famS:
                    geneUsedL[gene] = True

                # now set up familiesO to take this family and
                # determine the corresponding locusFamilies
                familiesO.initializeFamily(famNumCounter,familyMrca,[seedG1,seedG2])
                locusFamNumCounter,familiesO = createAllLocusFamiliesDescendingFromInternalNode(subtreeL,familyMrca,geneNames,famS,[seedG1,seedG2],famNumCounter,locusFamNumCounter,scoresO,minCoreSynThresh,minSynThresh,familiesO,familySubtreeNodeOrderL,nodeGenesL)
                famNumCounter+=1 # important to increment this after call to createAllLocusFamiliesDescendingFromInternalNode
                
    return geneUsedL,locusFamNumCounter,famNumCounter,familiesO

def createAllFamiliesAtTip(nodeGenesL,familyMrca,geneUsedL,tipFamilyRawThresholdD,scoresO,geneNames,minNormThresh,synAdjustThresh,synAdjustExtent,familiesO,minCoreSynThresh,minSynThresh,famNumCounter,locusFamNumCounter):
    '''Creates all Families and subsidiary LocusFamilies at the tip
given by familyMrca. Because we've come through the nodes in
pre-order, we know that all unused genes at this node are in a
families with mrca here. (they can still be multi gene families).
    '''

    unusedGenesAtThisTipS=set()
    for gene in nodeGenesL[familyMrca]: # familyMrca is a tip
        # gene is at this tip
        if not geneUsedL[gene]:
            # not used yet
            unusedGenesAtThisTipS.add(gene)

    # pull out the threshold we'll use given the strain
    tipFamilyRawThreshold = tipFamilyRawThresholdD[familyMrca]

    # Now we pull out families greedily. Simply take a gene,
    # and get all other genes that isSameFamily says are shared
    while len(unusedGenesAtThisTipS)>0:
        seed = unusedGenesAtThisTipS.pop()
        newFamS=set([seed])
        for newGene in unusedGenesAtThisTipS:

            if scoresO.isEdgePresentByEndNodes(seed,newGene):

                addIt = isSameFamily(seed,newGene,scoresO,geneNames,tipFamilyRawThreshold,minNormThresh,synAdjustThresh,synAdjustExtent)
                if addIt:
                    newFamS.add(newGene)

        # newFamS now contains a set of genes with significant
        # similarity. Remove it from unusedGenesAtThisTipS,
        # and create a new family from it.
        unusedGenesAtThisTipS.difference_update(newFamS)
        for gene in newFamS: # mark these as used
            geneUsedL[gene] = True


        familiesO.initializeFamily(famNumCounter,familyMrca)
        lfOL,locusFamNumCounter = createAllLocusFamiliesAtTip(newFamS,geneNames,familyMrca,scoresO,minCoreSynThresh,minSynThresh,famNumCounter,locusFamNumCounter)

        for lfO in lfOL:
            familiesO.addLocusFamily(lfO)

        famNumCounter+=1 # important to increment this after creating LocusFamilies

    return geneUsedL,locusFamNumCounter,famNumCounter,familiesO

def createLRSets(subtreeL,mrca,nodeGenesL,restrictS):
    '''At given mrca, obtain all genes in species in left branch and put
in leftS, and all genes from species in right branch to
rightS. Restrict each of these to be only genes in restrictS. If
restrictS is None, then use all genes.
    '''

    subtree=subtreeL[mrca]
    
    leftS=set()
    for tip in trees.leafList(subtree[1]):
        leftS.update(nodeGenesL[tip])
    rightS=set()
    for tip in trees.leafList(subtree[2]):
        rightS.update(nodeGenesL[tip])

    if restrictS != None:    
        leftS.intersection_update(restrictS)
        rightS.intersection_update(restrictS)
        
    return(leftS,rightS)

def closestMatch(gene,S,scoresO,geneNames,minNormThresh):
    '''Find the closest match to gene among the genes in the set S in the
graph scoresO. Eliminate any matches that have a norm score below
minNormThresh, a coreSynSc below minCoreSynThresh, or a synteny score
below minSynThresh.
    '''
    bestGene=None
    bestEdgeScore = -float('inf')
    for otherGene in scoresO.getConnectionsGene(gene):
        if otherGene in S:

            if isSameFamily(gene,otherGene,scoresO,geneNames,bestEdgeScore,minNormThresh,float('inf'),1):
                # the float('inf') is to avoid using synAdjustThresh here.
                bestEdgeScore = scoresO.getScoreByEndNodes(gene,otherGene,'rawSc')
                bestGene = otherGene
    return bestEdgeScore, gene, bestGene
    
def createSeedL(leftS,rightS,scoresO,geneNames,minNormThresh):
    '''Create a list which has the closest match for each gene on the
opposite side of the tree. e.g. if a gene is in tree[1] then we're
looking for the gene in tree[2] with the closest match. We eliminate
any matches that are below threshold for normalized or syntenty
scores. For each gene we get this closest match, put in a list, sort,
and return.
    '''
    seedL=[]
    for gene in leftS:
        seedL.append(closestMatch(gene,rightS,scoresO,geneNames,minNormThresh))
    for gene in rightS:
        seedL.append(closestMatch(gene,leftS,scoresO,geneNames,minNormThresh))
    seedL.sort(reverse=True)
    return seedL

def createFamilyFromSeed(g1,g2,geneUsedL,scoresO,leftS,rightS,geneNames,minRawThresh,minNormThresh,synAdjustThresh,synAdjustExtent):
    '''Based on a seed (seedScore, g1, g2) search for a family. Using the
PhiGs approach, we collect all genes which are closer to members of
the family than the two seeds are from each other. We have a normScore
threshold below which we will not add genes. We also have a synteny
adjustment of the raw score where we make the raw score between a pair
a bit better if their synteny is above synAdjustThresh. In general, if
a gene has syntenic connections to genes already in the family, this
makes us more confident that this gene belongs in the family. Returns
a set containing genes in the family.
    '''

    if geneUsedL[g1] or geneUsedL[g2]:
        # one of these has been used already, stop now.
        return None

    famS = set()
    genesToSearchForConnectionsS = set([g1,g2])

    while len(genesToSearchForConnectionsS) > 0:
        
        matchesS = getFamilyMatches(genesToSearchForConnectionsS,scoresO,leftS,rightS,famS,geneNames,minRawThresh,minNormThresh,synAdjustThresh,synAdjustExtent,geneUsedL)

        if matchesS == None:
            return None
        
        famS.update(genesToSearchForConnectionsS)
        genesToSearchForConnectionsS = matchesS
        
    return famS

def getFamilyMatches(genesToSearchForConnectionsS,scoresO,leftS,rightS,famS,geneNames,minRawThresh,minNormThresh,synAdjustThresh,synAdjustExtent,geneUsedL):
    ''''''
    matchesS=set()
    for famGene in genesToSearchForConnectionsS:
        for newGene in scoresO.getConnectionsGene(famGene):
            if newGene in leftS or newGene in rightS:
                # it is from a species descended from the node
                # we're working on
                if newGene not in genesToSearchForConnectionsS and newGene not in famS:
                    # it shouldn't been in our current list to search,
                    # or be one we've already put in the family (or
                    # we'll waste effort)
                    addIt = isSameFamily(famGene,newGene,scoresO,geneNames,minRawThresh,minNormThresh,synAdjustThresh,synAdjustExtent)
                    if addIt:
                        if geneUsedL[newGene]:
                            # this one's been used already. That
                            # means the whole family should be
                            # thrown out. Just stop now.
                            return None
                        else:
                            matchesS.add(newGene)
    return matchesS

def isSameFamily(famGene,newGene,scoresO,geneNames,minRawThresh,minNormThresh,synAdjustThresh,synAdjustExtent):
    '''Given famGene that is inside a family, and newGene we are
considering adding, check the various scores to determine if we should
add it. Return boolean.
    '''

    rawSc = scoresO.getScoreByEndNodes(famGene,newGene,'rawSc')
    normSc = scoresO.getScoreByEndNodes(famGene,newGene,'normSc')
    synSc = scoresO.getScoreByEndNodes(famGene,newGene,'synSc')
    
    addIt = False
    if normSc < minNormThresh:
        # doesn't meet norm score threshold for family
        # formation. (Modification of PhiGs)
        pass
    elif rawSc >= minRawThresh:
        # If its within the seed distance, add it
        # (basic PhiGs approach). we have the =
        # there in case minRawThresh is 1.
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
        if adjSc >= minRawThresh:
            addIt = True
                
    return addIt

def createAllLocusFamiliesDescendingFromInternalNode(subtreeL,familyMrca,geneNames,famGenesToSearchS,seedPairL,famNumCounter,locusFamNumCounter,scoresO,minCoreSynThresh,minSynThresh,familiesO,familySubtreeNodeOrderL,nodeGenesL):
    '''Given a family in famGenesToSearchS, break it up into subsidiary locus families
based on synteny. We iterate through the subtree rooted at familyMrca
in pre-order (ancestors first). Using seeds, we try to find groups
among famS that share high synteny.'''

    # split out LocusFamilies at non-syntenic locations
    for lfMrca,lchild,rchild in familySubtreeNodeOrderL:

        if lchild != None:
            # not a tip
            lfOL,locusFamNumCounter,famGenesToSearchS = createAllLocusFamiliesAtOneInternalNode(subtreeL,lfMrca,nodeGenesL,geneNames,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh,famNumCounter,locusFamNumCounter)
            
            for lfO in lfOL:
                familiesO.addLocusFamily(lfO)

        else:
            # we're at a tip.
            
            # get only the genes at this tip
            genesAtThisTipS = famGenesToSearchS.intersection(nodeGenesL[lfMrca])

            # remove them from famGenesToSearchS
            famGenesToSearchS.difference_update(genesAtThisTipS)

            # Get lf objects for all these genes
            lfOL,locusFamNumCounter = createAllLocusFamiliesAtTip(genesAtThisTipS,geneNames,lfMrca,scoresO,minCoreSynThresh,minSynThresh,famNumCounter,locusFamNumCounter)            

            # add to our families object
            for lfO in lfOL:
                familiesO.addLocusFamily(lfO)

    return locusFamNumCounter,familiesO

def createAllLocusFamiliesAtOneInternalNode(subtreeL,lfMrca,nodeGenesL,geneNames,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh,famNumCounter,locusFamNumCounter):
    '''Obtains all locus families at the internal node defined by lfMrca.'''

    lfOL = []
    while True:

        lfSeedPairL = createLFSeed(subtreeL,lfMrca,nodeGenesL,geneNames,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh)
        
        if lfSeedPairL == []:
            # there are no (more) seeds stradling this internal node,
            # break out
            break

        lfO,locusFamNumCounter,famGenesToSearchS = createLocusFamilyFromSeed(famNumCounter,locusFamNumCounter,lfMrca,lfSeedPairL,famGenesToSearchS,subtreeL,geneNames,scoresO,minCoreSynThresh,minSynThresh)
        lfOL.append(lfO)
        
    return lfOL,locusFamNumCounter,famGenesToSearchS
            
def createLFSeed(subtreeL,lfMrca,nodeGenesL,geneNames,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh):
    '''Given a set of genes famGenesToSearchS from a family, try to find a
seed based at lfMrca. A seed consists of two genes, one in the left
subtree and one in the right, which are syntenically consistent.
    '''
    
    leftS,rightS = createLRSets(subtreeL,lfMrca,nodeGenesL,famGenesToSearchS)
    
    for lGene in leftS:
        for rGene in rightS:
            if isSameLocusFamily(lGene,rGene,scoresO,geneNames,minCoreSynThresh,minSynThresh):
                return [lGene,rGene]
    return []
        
def createLocusFamilyFromSeed(famNumCounter,locusFamNumCounter,lfMrca,seedPairL,famGenesToSearchS,subtreeL,geneNames,scoresO,minCoreSynThresh,minSynThresh):
    '''Returns a LocusFamily object, containing genes associated with
those in seedPairL, in the subtree definied at lfMrca. Does single
linkage clustering, adding in anything in famGenesToSearchS with above
threshold synteny. Note that these seeds are not the seed from family formation (which might not be syntenic) but rather an independently generated pair which we know belong in the same LocusFamily
    '''

    lfO = LocusFamily(famNumCounter,locusFamNumCounter,lfMrca,[])
    locusFamNumCounter+=1
    
    for seedG in seedPairL:
        famGenesToSearchS.remove(seedG)
        lfO.addGene(seedG)

    subtree=subtreeL[lfMrca]
    strainL = trees.leafList(subtree)
        
    while True:
        genesToAddS = getLocusFamilyMatches(lfO,famGenesToSearchS,geneNames,strainL,scoresO,minCoreSynThresh,minSynThresh)
        if len(genesToAddS) == 0:
            break

        famGenesToSearchS.difference_update(genesToAddS)
        for gene in genesToAddS:
            lfO.addGene(gene)
        
    return lfO,locusFamNumCounter,famGenesToSearchS

def getLocusFamilyMatches(lfO,famGenesToSearchS,geneNames,strainL,scoresO,minCoreSynThresh,minSynThresh):
    '''Given a LocusFamily object lfO and some remaining genes, search
through the remaining genes to find those that match syntenically and
are in a child species of lfMrca. Return a list of genes to add.
    '''
    genesToAddS=set()
    for searchGene in famGenesToSearchS:
        # test if searchGene is in a child species of lfMrca
        if geneNames.numToStrainNum(searchGene) in strainL:
            # this searchGene is in a strain that is a child of the lfMrca we're working on
            for lfGene in lfO.genesL:
                # we don't use minNormThresh, minRawThresh, or
                # synAdjustThresh. If the pair have values above
                # minCoreSynThresh and minSynThres, then addIt will be
                # True.
                addIt = isSameLocusFamily(searchGene,lfGene,scoresO,geneNames,minCoreSynThresh,minSynThresh)
                if addIt:
                    genesToAddS.add(searchGene)
                    break
                
    return genesToAddS

def isSameLocusFamily(gene1,gene2,scoresO,geneNames,minCoreSynThresh,minSynThresh):
    '''Given two genes in the same family, determine if they meet the
synteny requirements to be put in the same LocusFamily. Returns
boolean.
    '''
    if not scoresO.isEdgePresentByEndNodes(gene1,gene2):
        # Within our families, there may be some gene-gene edges
        # missing due to the fact that blast could have just missed
        # significance etc. If the edge isn't there, then we do not
        # have evidence that these genes should be in the same locus
        # family, and we return false.
        
        # NOTE: An alternative approach would be to actually calculate
        # these scores here. But they're likely to be low...

        return False
    
    coreSynSc = scoresO.getScoreByEndNodes(gene1,gene2,'coreSynSc')
    
    if geneNames.isSameStrain(gene1,gene2):
        # same strain, only use core syneny
        if coreSynSc < minCoreSynThresh:
            # doesn't meet min core synteny requirement to be in same
            # LocusFamily
            addIt = False
        else:
            addIt = True
    else:
        # different strains, use both types of synteny score
        synSc = scoresO.getScoreByEndNodes(gene1,gene2,'synSc')
        if coreSynSc < minCoreSynThresh or synSc < minSynThresh:
            # one of the two types of synteny below threshold, so this
            # pair doesn't meet the requirements for being in the same
            # LocusFamily
            addIt = False
        else:
            addIt = True
            
    return addIt

def createAllLocusFamiliesAtTip(genesAtThisTipS,geneNames,lfMrca,scoresO,minCoreSynThresh,minSynThresh,famNumCounter,locusFamNumCounter):
    '''Given a set of genes famGenesToSearchS, search for all those found
on the tip given by lfMrca. Break these into LocusFamilies. Many will
be single gene LocusFamilies, but some may be multi-gene
    '''

    # Break these up into LocusFamilies. Many will be
    # single gene LocusFamilies, but some may be multi-gene
    lfGroupsL=[]
    while len(genesAtThisTipS) > 0:

        seed = genesAtThisTipS.pop()

        currentGroupS=set([seed])
        for gene in genesAtThisTipS:
            addIt = isSameLocusFamily(seed,gene,scoresO,geneNames,minCoreSynThresh,minSynThresh)
            if addIt:
                currentGroupS.add(gene)

        genesAtThisTipS.difference_update(currentGroupS)

        lfGroupsL.append(currentGroupS)

    lfOL=[]
    for lfGroupS in lfGroupsL:

        lfO = LocusFamily(famNumCounter,locusFamNumCounter,lfMrca,[])
        locusFamNumCounter+=1
        for gene in lfGroupS:
            lfO.addGene(gene)
        lfOL.append(lfO)

    return lfOL,locusFamNumCounter
    

def calcErrorScores(familyL,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD):
    '''NEEDS TO BE UPDATED. Given a list of family objects, calculate error scores for
each. These values get saved as attributes of the objects.
    '''
    for famO in familyL:
        famO.getPossibleErrorCt(scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD)


## Input/output

def writeFamilies(familiesO,geneNames,strainNum2StrD,fileName):
    '''Write all gene families to fileName, one family per line.'''
    f=open(fileName,'w')
    for fam in familiesO.iterFamilies():
        f.write(fam.fileStr(strainNum2StrD,geneNames)+'\n')
    f.close()


def readFamilies(familyFN,tree,geneNames,strainStr2NumD):
    '''Read the family file named familyFN, creating a Families object.
    '''

    familiesO = Families(tree)
    f=open(familyFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.split('\t')
        famNum=int(L[0])
        mrca = strainStr2NumD[L[1]]
        if L[2] == "-":
            seedPairL = None
        else:
            seedPairL = [L[2],L[3]]

        lfL = L[4:]

        familiesO.initializeFamily(famNum,mrca,seedPairL)

        for lfStr in lfL:
            lfSplitL = lfStr.rstrip().split(',')
            locusFamNum=int(lfSplitL[0])
            lfMrca = strainStr2NumD[lfSplitL[1]]
            geneL=[]
            for geneName in lfSplitL[2:]:
                geneL.append(geneNames.nameToNum(geneName))
            lfO = LocusFamily(famNum,locusFamNum,lfMrca,geneL)
            familiesO.addLocusFamily(lfO)

    f.close()
    return familiesO
