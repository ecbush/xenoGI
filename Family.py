
class Family:

    def __init__(self, idnum, mrca, genesL,numNodesInTree,geneNames,possibleErrorCt=None):
        '''Given an family ID number, a most recent common ancestor, and a
list of genes create a family object. Handles genes specified either
numerically or by name. numNodesInTree is the number of nodes found in
the tree (the full tree) we're working with for the
project. possibleErrorCt is a count of near miss genes, that either
were added but almost weren't, or were not added but almost were.
        '''

        self.id = idnum
        self.mrca = mrca

        # possibleErrorCt applies to multi-gene families. For single
        # gene families it is None.
        self.possibleErrorCt = possibleErrorCt 
        
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

        ##  TEMP MAKE NONE. FOR TESTING. CHANGE BACK LATER.
        #outL.append(str(self.possibleErrorCt))
        outL.append(str(None))


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


## USE BELOW (modified) for error count calculation

def isPossibleError(isInternal,g1,g2,scoresO,minNormThresh,minCoreSynThresh,minSynThresh,seedRawSc,synAdjustThresh,synAdjustExtent):
    '''Given g1 that is inside a family, and a g2 we are
considering adding, check the various scores to determine if we should
add it. isInternal is a boolean, telling us if the two genes are
within the same family or not. Return boolean.
    '''

    # get scores
    normSc = scoresO.getScoreByEndNodes(g1,g2,'normSc')
    synSc = scoresO.getScoreByEndNodes(g1,g2,'synSc')
    coreSynSc = scoresO.getScoreByEndNodes(g1,g2,'coreSynSc')


    # determine if it's a possible error and return magic numbers
    # here...but we'll be replacing this function with a more
    # systematic way run after family formation.
    if isPossibleError(isInternal,normSc,minNormThresh,0): # .5
        return 1
    elif isPossibleError(isInternal,coreSynSc,minCoreSynThresh,0): #.05
        return 1
    elif isPossibleError(isInternal,synSc,minSynThresh,0): # .5
        return 1
    else: return 0
    
def isPossibleErrorHelper(isInternal,sc,threshold,errorScoreDetermIncrement):
    '''Get better name. Determine if this was a near miss. If we added it, did we almost
    not do so? If we didn't add it, did we almost add it? If it's a
    near miss, it will contribute to our error score.
    '''
    if isInternal:
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

    
