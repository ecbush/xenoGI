
class Family:

    def __init__(self, idnum, mrca, genesL,numNodesInTree,geneNames):
        '''Given an family ID number, a most recent common ancestor, and a
list of genes create a family object. Handles genes specified either
numerically or by name. numNodesInTree is the number of nodes found in
the tree (the full tree) we're working with for the
project. possibleErrorCt is a count of near miss genes, that either
were added but almost weren't, or were not added but almost were.
        '''

        self.id = idnum
        self.mrca = mrca

        # create the family gene tuple, This has indexes corresponding
        # to nodes on the tree. At each position we have the genes at
        # that node.

        famGeneL = [[] for j in range(numNodesInTree)]
        geneNumL=[]
        for gene in genesL:
            # if we're being passed in genes with string names,
            # convert to number
            if type(gene)==int:
                geneNum=gene
            else:
                geneNum=geneNames.nameToNum(gene)

            strainNum=geneNames.numToStrainNum(geneNum)
            famGeneL[strainNum].append(geneNum)

        # tuple-ize
        newFamGeneL=[]
        for L in famGeneL:
            newFamGeneL.append(tuple(L))
        
        self.famGeneT = tuple(newFamGeneL)

    def getGeneNums(self):
        '''Extract and return all the gene numbers from this family.'''
        genesL=[]
        for geneT in self.famGeneT:
            genesL.extend(geneT)
        return genesL

    def getGeneNames(self,geneNames):
        '''Extract and return all the gene names from this family.'''
        genesL=[]
        for num in self.getGeneNums():
            genesL.append(geneNames.numToName(num))
        return genesL

    def isInStrain(self,strainNum):
        '''Is this family present in strainNum, return boolean.'''
        if len(self.famGeneT[strainNum]) > 0:
            return True
        else: return False
        
    def __repr__(self):
        '''String representation of a family containing family number.'''
        return "<Family: "+str(self.id) + ">"

    def fileStr(self,geneNames,strainNum2StrD):
        '''String representation suitable for writing the family to
file. Genes and mrca are expressed in word form.'''

        outL =[str(self.id)]
        outL.append(strainNum2StrD[self.mrca])
        outL.extend(self.getGeneNames(geneNames))
        return "\t".join(outL)

    def getOutsideConnections(self,scoresO):
        '''Given a score object, return a tuple of all outside genes with
connections to this family.'''
        allGenesInFamL = self.getGeneNums()
        otherGenesS=set()
        for geneNum in allGenesInFamL:
            for otherGene in scoresO.getConnectionsGene(geneNum):
                if not otherGene in allGenesInFamL:
                    otherGenesS.add(otherGene)
        return tuple(otherGenesS)
    
    def getPossibleErrorCt(self, scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD):
        '''Given a family, returns the number of near misses internally and externally'''
        internalPossibleErrors,externalPossibleErrors = 0,0
        allGenesInFamL = self.getGeneNums()
        internalEdgeL = makeMSN(allGenesInFamL,scoresO)
        externalGenesTuple=self.getOutsideConnections(scoresO)
        externalEdgeL = makeExternalEdgeL(externalGenesTuple, allGenesInFamL, scoresO)
        for g1,g2 in internalEdgeL:
            internalPossibleErrors+=isPossibleErrorInternal(g1,g2,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD)
        for g1,g2 in externalEdgeL:
            externalPossibleErrors+=isPossibleErrorExternal(g1,g2,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD)

        self.possibleErrorCt = internalPossibleErrors + externalPossibleErrors


# Functions to help in calculating possible error count

def makeExternalEdgeL(externalGenesTuple, genesInFam, scoresO):
    '''given a family, returns a list of tuples '''
    returnEdgeL = [0]*len(externalGenesTuple)
    for extGeneIndex in range(0,len(externalGenesTuple)):
        bestScore = float('-inf')
        extGene = externalGenesTuple[extGeneIndex]
        connectedGenes = scoresO.getConnectionsGene(extGene)
        for conGene in connectedGenes:
            if conGene in genesInFam:
                rawSc=scoresO.getScoreByEndNodes(extGene,conGene, 'rawSc')
                if rawSc>bestScore:
                    bestScore=rawSc
                    returnEdgeL[extGeneIndex]=(extGene,conGene)
    return returnEdgeL

def makeInternalEdgeL(genesInFam,scoresO):
    '''Given a list of all genes in the family, it returns a
list of lists sorted in descending order by raw score'''
    # compare every gene to every other gene that comes
    # in the list after it
    # format: [[rawScore, gene1, gene2], ....]] for every edge in the fam
    totalGenesNum = len(genesInFam)
    internalEdgeL = []
    for gene1Index in range(0,totalGenesNum-1):
        for gene2Index in range(gene1Index+1,totalGenesNum):
            g1=genesInFam[gene1Index]
            g2=genesInFam[gene2Index]
            if scoresO.isEdgePresentByEndNodes(g1,g2):
                rawSc = scoresO.getScoreByEndNodes(g1,g2, 'rawSc')
                internalEdgeL.append([rawSc,g1,g2])
    internalEdgeL.sort(reverse=True)
    return internalEdgeL


def makeUnconnectedGsetL(genesInFam):
    '''Given a list of all genes in the family, make a set of disjoint genes'''
    unconnectedGsetL = []
    for gene in genesInFam:
        geneSet = set()
        geneSet.add(gene)
        unconnectedGsetL.append(geneSet)
    return unconnectedGsetL


def updateConnections(gSetL, g1, g2):
    '''Given 2 genes and a gene set list, returns True if the genes
are disconnected, false otherwise. If genes are disconnected, join the
genes and update gSetL accordingly'''
    for setIndex1 in range(0,len(gSetL)):
        if g1 in gSetL[setIndex1]:
               if g2 in gSetL[setIndex1]: return False,gSetL
               else:
                   for setIndex2 in range(0, len(gSetL)):
                       if g2 in gSetL[setIndex2]:
                           gSetL[setIndex1]=gSetL[setIndex1].union(gSetL[setIndex2])
                           gSetL.remove(gSetL[setIndex2])
                           return True, gSetL
    return False, gSetL
        

def makeMSN(genesInFam,scoresO):
    '''Given a tuple of genes in the family, it makes the max spanning network based on
the largest raw scores'''
    returnEdgeList = []
    internalEdgeL = makeInternalEdgeL(genesInFam,scoresO)
    gSetL = makeUnconnectedGsetL(genesInFam)
    while len(gSetL) > 1:
        currEdge = internalEdgeL[0]
        internalEdgeL.remove(currEdge)
        g1,g2 = currEdge[1],currEdge[2]
        isDisconnected, gSetL = updateConnections(gSetL, g1, g2)
        if isDisconnected:
            returnEdgeList.append((g1,g2))
    return returnEdgeList

def isPossibleErrorInternal(g1,g2,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD):
    '''Given g1 and g2 both inside a family, check the various scores to
    determine if it was a near miss, ie almost was not in the
    family. Returns boolean.
    '''
    # get scores
    normSc = scoresO.getScoreByEndNodes(g1,g2,'normSc')
    synSc = scoresO.getScoreByEndNodes(g1,g2,'synSc')
    coreSynSc = scoresO.getScoreByEndNodes(g1,g2,'coreSynSc')
    
    if not (normSc - famErrorScoreIncrementD['normSc']) >= minNormThresh:
        # Was over, but now is below threshold. Is a possible error
        return True
    elif not (coreSynSc - famErrorScoreIncrementD['coreSynSc']) >= minCoreSynThresh:
        return True
    elif not (synSc - famErrorScoreIncrementD['synSc']) >= minSynThresh:
        return True
    else: return False

def isPossibleErrorExternal(g1,g2,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD):
    '''Given g1 inside and g2 outside a family, check the various scores
    to determine if it was a near miss, ie almost was put in the
    family. We consider normSc, synSc, and coreSynSc. We define near
    miss to be cases when 2/3 of these were above threshold, and the
    third was below, but within an increment of the threshold
    (increments provided in famErrorScoreIncrementD). Returns boolean.
    '''
    # get scores
    normSc = scoresO.getScoreByEndNodes(g1,g2,'normSc')
    synSc = scoresO.getScoreByEndNodes(g1,g2,'synSc')
    coreSynSc = scoresO.getScoreByEndNodes(g1,g2,'coreSynSc')

    # w/out increment
    normB = normSc >= minNormThresh
    coreSynB = coreSynSc >= minCoreSynThresh
    synB = synSc >= minSynThresh

    # with increment
    normIncB = (normSc + famErrorScoreIncrementD['normSc']) >= minNormThresh
    coreSynIncB = (coreSynSc + famErrorScoreIncrementD['coreSynSc']) >= minCoreSynThresh
    synIncB = (synSc + famErrorScoreIncrementD['synSc']) >= minSynThresh
    
    if coreSynB and synB and not normB and normIncB:
        # core and syn were over without increment, but norm is only
        # over once we add the increment
        return True
    elif normB and synB and not coreSynB and coreSynIncB:
        return True
    elif coreSynB and normB and not synB and synIncB:
        return True
    else: return False
