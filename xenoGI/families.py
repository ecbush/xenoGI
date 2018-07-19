# Functions for a modiefied version of the PhiGs algorithm
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# (We've added the use of synteny)
import sys
from . import trees
from .Family import *


## Main function



def createFamiliesO(tree,scoresO,geneNames,subtreeL,minNormThresh,minCoreSynThresh,minSynThresh,synAdjustThresh,synAdjustExtent,outputSummaryF,strainNum2StrD,familyFN):
    '''Given a graph of genes and their similarity scores find families
using a PhiGs-like algorithm, with synteny also considered.
    '''

    nodeOrderL=createNodeProcessOrderList(tree)

    # initialize scoresO.nodeConnectL for use below
    scoresO.createNodeConnectL(geneNames)
    
    geneUsedL = [False for x in geneNames.nums]

    # create an object of class Families to store this in.
    familiesO = Families(tree)
    famNumCounter = 0
    locusFamNumCounter = 0
    
    for familyMrca,lchild,rchild in nodeOrderL:
        if lchild != None:
            # not a tip
            subtree=subtreeL[familyMrca]

            # for use in createSubFamilies call below
            familySubtreeNodeOrderL = createNodeProcessOrderList(subtree)
            
            leftS,rightS,outgroupS = createLRSets(subtree,geneNames,geneNames.nums)
            seedL = createSeedL(leftS,rightS,scoresO,minNormThresh,minCoreSynThresh,minSynThresh)
            for seed in seedL:
                # each seed corresponds to a prospective gene family.
                seedRawSc,seedG1,seedG2 = seed

                if seedRawSc == -float('inf'):
                    # we've gotten to the point in the seed list with
                    # genes having no match on the other branch
                    break
                else:
                    # getting initial family, using only raw score and synteny bump
                    famS=getFamily(seedRawSc,seedG1,seedG2,scoresO,leftS,rightS,minNormThresh,synAdjustThresh,synAdjustExtent,geneUsedL)
                    
                    if famS == None:
                        # one of the genes the the family was already
                        # used, so getFamily returned None
                        continue
                    else:
                        # none of the genes in famS used yet
                        for gene in famS:
                            geneUsedL[gene] = True

                        # now set up familiesO to take this family and
                        # determine the corresponding locusFamilies
                        
                        familiesO.initializeFamily(famNumCounter,familyMrca,[seedG1,seedG2])
                        famNumCounter,locusFamNumCounter,familiesO = createLocusFamilies(subtreeL,familyMrca,geneNames,famS,[seedG1,seedG2],famNumCounter,locusFamNumCounter,scoresO,minCoreSynThresh,minSynThresh,familiesO,familySubtreeNodeOrderL)
                        famNumCounter+=1



    genesInMultiGeneFamsS = familiesO.getAllGenes()
    

    multiGeneFamNum=famNumCounter
    print("Number of multigene families",multiGeneFamNum,file=outputSummaryF)

    # Add on remaining genes as single gene families.
    for gene in geneNames.nums:
        if not gene in genesInMultiGeneFamsS:
            strain = geneNames.numToStrainNum(gene)
            lfO = LocusFamily(famNumCounter,locusFamNumCounter,strain)
            familiesO.initializeFamily(famNumCounter,strain)
            familiesO.addLocusFamily(lfO)
            locusFamNumCounter+=1
            famNumCounter+=1
            

    print("Number of single gene families",famNumCounter-multiGeneFamNum,file=outputSummaryF)
    print("Number of total families",famNumCounter,file=outputSummaryF)

    writeFamilies(familyL,geneNames,strainNum2StrD,familyFN)

    return familiesO

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


def createLRSets(tree,geneNames,geneListL):
    '''For every gene in our data, put it into one of three sets. Left,
right, or outgroup. Genes in the left set are found in a species on
the left branch of tree.'''

    leftSpeciesS=set(trees.leafList(tree[1]))
    rightSpeciesS=set(trees.leafList(tree[2]))
    
    leftS=set()
    rightS=set()
    outgroupS=set()
    for geneNum in geneListL: # all genes
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

            if isFamily(gene,otherGene,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,bestEdgeScore,float('inf'),1):
                # the float('inf') is to avoid using synAdjustThresh here.
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

def getFamily(seedRawSc,g1,g2,scoresO,leftS,rightS,minNormThresh,synAdjustThresh,synAdjustExtent,geneUsedL):
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

    ## Maybe re-write this with a helper function, analagous to
    ## getLocusFamilyMatches? That approach is clearer if it works
    ## here.
    
    if geneUsedL[g1] or geneUsedL[g2]:
        # one of these has been used already, stop now.
        return None

    alreadySearchedS = set()
    notYetSearchedS = set([g1,g2])

    while len(notYetSearchedS) > 0:
        newGenesS=set()
        for famGene in notYetSearchedS:
            for newGene in scoresO.getConnectionsGene(famGene):
                if newGene in leftS or newGene in rightS:
                    # it is from a species descended from the node
                    # we're working on
                    if newGene not in alreadySearchedS and newGene not in newGenesS and newGene not in notYetSearchedS:
                        addIt = isFamily(famGene,newGene,scoresO,minNormThresh,-float('inf'),-float('inf'),seedRawSc,synAdjustThresh,synAdjustExtent)
                        if addIt:
                            if geneUsedL[newGene]:
                                # this one's been used already. That
                                # means the whole family should be
                                # thrown out. Just stop now.
                                return None
                            else:
                                newGenesS.add(newGene)
                    
            alreadySearchedS.add(famGene)
        notYetSearchedS = newGenesS
    return alreadySearchedS


def isFamily(famGene,newGene,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,minRawThresh,synAdjustThresh,synAdjustExtent):
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

def createLocusFamilies(subtreeL,familyMrca,geneNames,famS,seedPairL,famNumCounter,locusFamNumCounter,scoresO,minCoreSynThresh,minSynThresh,familiesO,familySubtreeNodeOrderL):
    '''Given a family in famS, break it up into subsidiary locus families
based on synteny. We iterate through the subtree rooted at familyMrca
in pre-order (ancestors first). Using seeds, we try to find groups
among famS that share high synteny.'''

    #print("seedPairL",seedPairL)

    subtree=subtreeL[familyMrca]
    leftS,rightS,outgroupS = createLRSets(subtree,geneNames,famS)

    famGenesToSearchS = leftS.copy()
    famGenesToSearchS.update(rightS)

    for seedG in seedPairL:
        famGenesToSearchS.remove(seedG)

    # get the LocusFamily associated with the original seed
    lfO,locusFamNumCounter,famGenesToSearchS = createLocusFamilyFromSeed(famNumCounter,locusFamNumCounter,familyMrca,seedPairL,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh)

    familiesO.addLocusFamily(lfO)
    
    # get any remaining LocusFamilies at non-syntenic locations
    
    for lfMrca,lchild,rchild in familySubtreeNodeOrderL:

        if lchild != None:
            # not a tip
            subtree=subtreeL[lfMrca]
            leftS,rightS,outgroupS = createLRSets(subtree,geneNames,famGenesToSearchS)

            # not using minNormThresh
            seedL = createSeedL(leftS,rightS,scoresO,-float('inf'),minCoreSynThresh,minSynThresh)

            for seed in seedL:
                seedRawSc,g1,g2 = seed

                if seedRawSc == -float('inf'):
                    # we've gotten to the point in the seed list with
                    # genes having no match on the other branch
                    break
                else:
                    startGeneL = [g1,g2]
                    lfO,locusFamNumCounter,famGenesToSearchS = createLocusFamilyFromSeed(famNumCounter,locusFamNumCounter,lfMrca,startGeneL,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh)
                    familiesO.addLocusFamily(lfO)

    # any remaining genes must be on tips, and belong in their own locusFamily
    for gene in famGenesToSearchS:
        strain=geneNames.numToStrainNum(gene)
        lfO = LocusFamily(famNumCounter,locusFamNumCounter,strain)
        locusFamNumCounter+=1
        lfO.add(gene)
        familiesO.addLocusFamily(lfO)
        
    

    return famNumCounter,locusFamNumCounter,familiesO
    

def createLocusFamilyFromSeed(famNumCounter,locusFamNumCounter,lfMrca,startGeneL,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh):
    '''Returns a LocusFamily object, containing genes associated with
those in startGeneL, in the subtree definied at lfMrca
(famGenesToSearchS is already restricted to the subtree defined by
lfMrca). Does single linkage clustering, adding in anything in
famGenesToSearchS with above threshold synteny.
    '''
    
    lfO = LocusFamily(famNumCounter,locusFamNumCounter,lfMrca)
    locusFamNumCounter+=1

    for gene in startGeneL:
        lfO.addGene(gene)

    while True:
        genesToAddS = getLocusFamilyMatches(lfO,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh)
        if len(genesToAddS) == 0:
            break

        famGenesToSearchS.difference_update(genesToAddS)
    
    return lfO,locusFamNumCounter,famGenesToSearchS

def getLocusFamilyMatches(lfO,famGenesToSearchS,scoresO,minCoreSynThresh,minSynThresh):
    '''Given a LocusFamily object lfO and some remaining genes, search through the remaining genes to find those that match syntenically. Return a list of genes to add.'''
    
    genesToAddS=set()
    
    for lfGene in lfO.genesL:
        for searchGene in famGenesToSearchS:
            # we don't use minNormThresh, minRawThresh, or
            # synAdjustThresh. If the pair have values above
            # minCoreSynThresh and minSynThres, then addIt will be
            # True.
            addIt = isFamily(lfGene,searchGene,scoresO,-float('inf'),minCoreSynThresh,minSynThresh,-float('inf'),float('inf'),1)
            if addIt:
                genesToAddS.add(searchGene)
                
    return genesToAddS
                             

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
    for fam in familiesO.iterFamily():
        f.write(fam.fileStr(geneNames,strainNum2StrD)+'\n')
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
            lfSplitL = lfStr.split(',')
            locusFamNum=int(lfSplitL[0])
            lfMrca = strainStr2NumD[lfSplitL[1]]
            geneL=[]
            for geneName in lfSplitL[2:]:
                geneL.append(geneNames.nameToNum(geneName))
            lfO = LocusFamily(locusFamNum,famNum,lfMrca)
            familiesO.addLocusFamily(lfO)

    f.close()
    return familiesO
