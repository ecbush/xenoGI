# Functions for a modiefied version of the PhiGs algorithm
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# (We've added the use of synteny)
# Eliot Bush 8/2016
import sys
import trees
from Family import *

def nodeDetails(tree,timeOnBrLeadingHere):
    '''Given a tree, return list of (time,node #, left node #, right node
#) for each node. timeOnBrLeadingHere gives the time which should be
added to the time values in the tuple we're making.
    '''
    if tree[1] == ():
        return [(tree[3] + timeOnBrLeadingHere , tree[0], None, None)]
    else:
        l = nodeDetails(tree[1], tree[3] + timeOnBrLeadingHere)
        r = nodeDetails(tree[2], tree[3] + timeOnBrLeadingHere)
        return l + r + [(tree[3] + timeOnBrLeadingHere , tree[0], tree[1][0], tree[2][0])]

        
def createNodeProcessOrderList(tree):
    '''Given a tree, output a list specifying the order of nodes we should
process in the PhiGs algorithm. For each node we give a tuple
(time,node #, left node #, right node #)
    '''
    L=nodeDetails(tree,0)
    L.sort()
    return L

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

def closestMatch(gene,S,simG,synScoresG,minSynThresh):
    '''Find the closest match to gene among the genes in the set S in the
graph simG. Eliminate any matches that have a synteny score below
minSynThresh
    '''
    bestGene=None
    bestEdgeScore = -float('inf')
    for edge in simG.edges_iter(gene):
        if edge[1] in S:
            if simG.get_edge_data(*edge)['score'] > bestEdgeScore:
                if synScoresG.get_edge_data(*edge)['score'] > minSynThresh:
                    bestEdgeScore =  simG.get_edge_data(*edge)['score']
                    bestGene = edge[1]
    return bestEdgeScore, gene, bestGene
    
def createSeedL(leftS,rightS,simG,synScoresG,minSynThresh):
    '''Create a list which has the closest match for each gene on the
opposite side of the tree. e.g. if a gene is in tree[1] then we're
looking for the gene in tree[2] with the closest match. We eliminate
any matches that have a synteny score below minSynThresh. For each
gene we get this closest match, put in a list, sort, and return.
    '''
    seedL=[]
    for gene in leftS:
        seedL.append(closestMatch(gene,rightS,simG,synScoresG,minSynThresh))
    for gene in rightS:
        seedL.append(closestMatch(gene,leftS,simG,synScoresG,minSynThresh))
    seedL.sort(reverse=True)
    return seedL

def getFamily(seedSimScore,g1,g2,simG,synScoresG,leftS,rightS,minSynThresh,synAdjustThresh,synAdjustMaxExtent):
    '''Based on a seed (seedScore, g1, g2) search for a family. Using the
PhiGs approach, we collect all genes which are closer to members of
the family than the two seeds are from each other. But, we also use
synteny, by using synteny scores to adjust the sim scores between
genes we are considering adding. In general, if a gene has syntenic
connectsions to genes already in the family, this makes us more
confident that this gene belongs in the family. Returns a set
containing genes in the family.
    '''
    alreadySearchedS = set()
    notYetSearchedS = set([g1,g2])

    while len(notYetSearchedS) > 0:
        newGenesS=set()
        for gene in notYetSearchedS:
            for edge in simG.edges_iter(gene):
                newGene=edge[1]
                if newGene in leftS or newGene in rightS:
                    # it is from a species descended from the node we're working on
                    if newGene not in alreadySearchedS and newGene not in newGenesS:
                        sc = simG.get_edge_data(*edge)['score']
                        synsc = synScoresG.get_edge_data(*edge)['score']
                        if synsc < minSynThresh:
                            # doesn't meet min synteny requirement for family formation
                            continue
                        elif sc >= seedSimScore:
                            # if its above the regular score
                            # threshold, add it. we have the = there
                            # in case seedSimScore is 1.
                            newGenesS.add(newGene)
                        elif synsc > synAdjustThresh:
                            # its above the syn score adjustment
                            # threshold, use it to adjust the
                            # score. Adjustment magnitude depends linearly
                            # on synsc, and ranges between 0 and
                            # synAdjustMaxExtent.
                            adjustment = synAdjustMaxExtent * ( synsc - synAdjustThresh ) / (1.0 - synAdjustThresh)
                            adjsc = sc + adjustment
                            if adjsc > 1: adjsc = 1 # truncate back to 1
                            if adjsc > seedSimScore:
                                newGenesS.add(newGene)
                        else: continue # nothing worked, don't add it.
                            
            alreadySearchedS.add(gene)
        notYetSearchedS = newGenesS
    return alreadySearchedS
                
    
def families(tree,subtreeL,geneNames,simG,synScoresG,minSynThresh,synAdjustThresh,synAdjustMaxExtent,familyFN,strainNum2StrD):
    '''Given a graph of genes and their similarity scores (simG) find
families using a PhiGs-like algorithm, with synteny also considered.'''

    nodeOrderL=createNodeProcessOrderList(tree)
    
    geneUsedL = [False for x in geneNames.nums]
    
    multiGeneFamilyL=[]
    
    for t,node,lnode,rnode in nodeOrderL:
        if lnode != None:
            # not a tip
            subtree=subtreeL[node]
            leftS,rightS,outgroupS = createLRSets(subtree,geneNames)
            seedL = createSeedL(leftS,rightS,simG,synScoresG,minSynThresh)
            #print(len(seedL),len(leftS),len(rightS),len(outgroupS),'xx')
            for seed in seedL:
                seedSimScore,g1,g2 = seed
                #print(seedSimScore,g1,g2,'xx')

                if seedSimScore == -float('inf'):
                    # we've gotten to the point in the seed list with
                    # genes having no match on the other branch
                    break
                elif synScoresG.get_edge_data(g1,g2)['score'] < minSynThresh:
                    # doesn't meet min synteny requirement for family formation
                    continue
                else:
                    family=getFamily(seedSimScore,g1,g2,simG,synScoresG,leftS,rightS,minSynThresh,synAdjustThresh,synAdjustMaxExtent)
                    
                    if any((geneUsedL[gene] for gene in family)):
                        continue
                    else:
                        # none of the genes in this family used yet
                        for gene in family:
                            geneUsedL[gene] = True
                        multiGeneFamilyL.append((node,family))

    # Create Family objects 
    genesInMultiGeneFamsS = set()
    familyL=[]
    famNum = 0
    for mrca,famS in multiGeneFamilyL:
        genesInMultiGeneFamsS.update(famS) # add all genes in fam
        familyL.append(Family(famNum,mrca,famS,trees.nodeCount(tree),geneNames))
        famNum+=1

    multiGeneFamNum=famNum
    print("Number of multigene families",multiGeneFamNum,file=sys.stderr)

    # Add on remaining genes as single gene families.
    for gene in geneNames.nums: 
        if not gene in genesInMultiGeneFamsS:
            mrca = geneNames.numToStrainNum(gene)
            familyL.append(Family(famNum,mrca,[gene],trees.nodeCount(tree),geneNames))
            famNum+=1

    print("Number of single gene families",famNum-multiGeneFamNum,file=sys.stderr)
    print("Number of total families",famNum,file=sys.stderr)
    
    writeFamilies(familyL,geneNames,strainNum2StrD,familyFN)

    return tuple(familyL)


def writeFamilies(familyL,geneNames,strainNum2StrD,fileName):
    '''Write all gene families to fileName, one family per line.'''
    f=open(fileName,'w')
    for fam in familyL:
        f.write(fam.fileStr(geneNames,strainNum2StrD)+'\n')
    f.close()


def readFamilies(familyFN,tree,geneNames,strainStr2NumD):
    '''Read the family file named familyFN, creating a list of family
objects.'''
    
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
        genesL = L[2:]

        famL.append(Family(famNum,mrca,genesL,trees.nodeCount(tree),geneNames))
    return tuple(famL)
