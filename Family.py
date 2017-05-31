
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
        # to nodes on the tree. At each position we have another tuple
        # (gene count, (tuple of genes))

        familyL = [[0,[]] for j in range(numNodesInTree)]
        geneNumL=[]
        for gene in genesL:
            # if we're being passed in genes with string names,
            # convert to number
            if type(gene)==int:
                geneNum=gene
            else:
                geneNum=geneNames.nameToNum(gene)

            geneNumL.append(geneNum)
            strainNum=geneNames.numToStrainNum(geneNum)
            familyL[strainNum][0]+=1
            familyL[strainNum][1].append(geneNum)

        # tuple-ize
        newFamilyL=[]
        for ct,L in familyL:
            newFamilyL.append((ct,tuple(L)))
        
        self.famGeneT = tuple(newFamilyL)
        self.geneT = tuple(geneNumL)
        
    def __repr__(self):
        '''String representation of a family containing family number.'''
        return "<Family: "+str(self.id) + ">"

    def fileStr(self,geneNames,strainNum2StrD):
        '''String representation suitable for writing the family to
file. Genes and mrca are expressed in word form.'''

        outL =[str(self.id)]
        outL.append(strainNum2StrD[self.mrca])

        for ct,geneT in self.famGeneT:
            for gene in geneT:
                outL.append(geneNames.numToName(gene))
        return "\t".join(outL)

    def getOutsideConnections(self,scoresO):
        '''Given a score object, return a tuple of all outside genes with
connections to this family.'''
        otherGenesS=set()
        for geneNum in self.geneT:
            for otherGene in scoresO.getConnectionsGene(geneNum):
                if not otherGene in self.geneT:
                    otherGenesS.add(otherGene)
        return tuple(otherGenesS)
    
    def getPossibleErrorCt(self, scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD):
        '''Given a family, returns the number of near misses internally and externally'''
        internalPossibleErrors,externalPossibleErrors = 0,0
        allGenesInFam = self.geneT
        internalEdgeL = makeMSN(allGenesInFam,scoresO)
        externalGenesTuple=self.getOutsideConnections(scoresO)
        externalEdgeL = makeExternalEdgeL(externalGenesTuple, allGenesInFam, scoresO)
        for g1,g2 in internalEdgeL:
            internalPossibleErrors+=isPossibleErrorInternal(g1,g2,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD)
        for g1,g2 in externalEdgeL:
            externalPossibleErrors+=isPossibleErrorExternal(g1,g2,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,famErrorScoreIncrementD)

        self.internalPossibleErrorCt = internalPossibleErrors
        self.externalPossibleErrorCt = externalPossibleErrors


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
    family. Returns boolean.
    '''
    # get scores
    normSc = scoresO.getScoreByEndNodes(g1,g2,'normSc')
    synSc = scoresO.getScoreByEndNodes(g1,g2,'synSc')
    coreSynSc = scoresO.getScoreByEndNodes(g1,g2,'coreSynSc')

    
    if (normSc + famErrorScoreIncrementD['normSc']) >= minNormThresh:
        # Was below, but now is over threshold. Is a possible error
        return True
    elif (coreSynSc + famErrorScoreIncrementD['coreSynSc']) >= minCoreSynThresh:
        return True
    elif (synSc + famErrorScoreIncrementD['synSc']) >= minSynThresh:
        return True
    else: return False
