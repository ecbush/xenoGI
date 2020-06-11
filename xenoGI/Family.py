
class LocusFamily:
    def __init__(self, famNum, locusFamNum, lfMrca,locusNum=None):
        '''Initialize a LocusFamily object with a family number, a LocusFamily
number, and the mrca for this locus family (which may differ from the
mrca for its family.'''
        self.famNum = famNum
        self.locusFamNum = locusFamNum # a unique number for each lf
        self.lfMrca = lfMrca
        self.locusNum = locusNum # indicates location, several lf may share
        self.geneD = {}

    def addGene(self, gene,genesO):
        '''Add a gene (in numerical form) to this LocusFamily.'''
        strain = genesO.numToStrainName(gene)
        if strain in self.geneD:
            self.geneD[strain].append(gene)
        else:
            self.geneD[strain] = [gene]
        
    def addGenes(self,genesL,genesO):
        '''Add a list (or set etc.) of genes to this LocusFamily.'''
        for gene in genesL:
            self.addGene(gene,genesO)

    def iterGenes(self):
        '''Iterate through all genes in this locus family.'''
        for L in self.geneD.values():
            for gene in L:
                yield gene

    def iterGenesByStrain(self,strain):
        '''Iterate through all genes in this locus family that are in a
        particular strain.
        '''
        if strain in self.geneD:
            for gene in self.geneD[strain]:
                yield gene
        return

    def iterStrains(self):
        '''Iterator that yields all the strains in which this locus family
occurs.'''
        for strain in self.geneD.keys():
            yield strain

    def origin(self,familiesO):
        '''Determine the event or events at the origin of this locus
family.
        '''
        fam = familiesO.getFamily(self.famNum)
        return fam.reconD[self.lfMrca]
        
    def getStr(self,genesO,sep):
        '''Return a string representation of a single LocusFamily. Separator
between elements given by sep. Elements are: locusFamNum lfMrca gene1
gene2...
        '''
        outL=[str(self.locusFamNum),self.lfMrca,str(self.locusNum)]
        
        for geneNum in self.iterGenes():
            outL.append(genesO.numToName(geneNum))

        return sep.join(outL)
    
    def fileStr(self,genesO):
        '''Return a string representation of a single LocusFamily. Format is
        comma separated: 
        '''
        return self.getStr(genesO,',')

    def __repr__(self):
        '''String representation of a LocusFamily, for display purposes.'''
        return "<lf "+str(self.locusFamNum)+">"
        
        
class Family:
    def __init__(self,famNum,mrca,geneTree=None,reconD=None,sourceFam=None):
        '''Initialize an object of class Family.'''

        self.famNum = famNum
        self.mrca = mrca
        self.locusFamiliesL = []   # will contain locusFamily objects for this family
        self.geneTree = geneTree   # rooted gene tree object
        self.reconD = reconD       # dict giving mapping of gene tree onto species tree
        self.sourceFam = sourceFam # only for originFams. gives ifam num that originFam came from
        
    def addLocusFamily(self,lfO):
        self.locusFamiliesL.append(lfO)

    def getLocusFamilies(self):
        return self.locusFamiliesL

    def iterGenes(self):
        '''Iterate through all the genes associated with this family in numerical form.'''
        for lfO in self.getLocusFamilies():
            for gene in lfO.iterGenes():
                yield gene

    def iterStrains(self):
        '''Iterator that yields all the strains in which this family
occurs.
        '''
        for lfO in self.getLocusFamilies():
            for strain in lfO.iterStrains():
                yield strain
                
    def getOutsideConnections(self,scoresO):
        '''Given a score object, return a set of all outside genes with
connections to this family.'''
        allGenesInFamS = set(self.iterGenes())
        otherGenesS=set()
        for geneNum in allGenesInFamS:
            for otherGene in scoresO.getConnectionsGene(geneNum):
                if not otherGene in allGenesInFamS:
                    otherGenesS.add(otherGene)
        return otherGenesS

    def addGeneTree(self,geneTree):
        '''Given tree object geneTree, store as attribute. Currently this
object is assumed to be a tuple tree.'''
        self.geneTree = geneTree

    def addReconciliation(self,rD):
        '''Given a dictionary rD storing reconciliation values, store as attribute.'''
        self.reconD = rD
        
    def fileStr(self,genesO):
        '''Return string representation of single family. Format is: famNum <tab> 
        mrca <tab> locusFamNum1,locusFamGenes <tab>
        locusFamNum2,locusFamGenes... <tab> geneTree <tab> reconciliation
        The LocusFamily object representations are comma separated.
        The geneTree is string version of a tuple.
        reconciliation is a string of the dict.
        In future improve this...
        '''

        # family number and mrca
        outL =[str(self.famNum),self.mrca]

        # geneTree
        if self.geneTree != None:
            outL.append(str(self.geneTree))
        else:
            outL.append("None")

        # reconciliation
        if self.reconD != None:
            outL.append(str(self.reconD))
        else:
            outL.append("None")

        # source family
        if self.sourceFam != None:
            outL.append(str(self.sourceFam))
        else:
            outL.append("None")
            
        # locus families
        for lfO in self.locusFamiliesL:
            outL.append(lfO.fileStr(genesO))

        return "\t".join(outL)

    def __repr__(self):
        return "<fam:"+str(self.famNum)+">"
    

class Families:

    def __init__(self,speciesTree):
        '''Initialize an object of class Families.'''

        self.speciesTree = speciesTree
        self.locusFamiliesD = {}
        self.familiesD = {} 

        ## locusFamiliesD has key locusFamNum and value a LocusFamily
        ## object. familiesD has key famNum and value
        ## a Family object.
        
    def initializeFamily(self,famNum,mrca,geneTree=None,reconD=None,sourceFam=None):
        '''Set up an entry for family famNum.'''
        # the seed genes are the original PHiGs seed.

        self.familiesD[famNum] = Family(famNum,mrca,geneTree,reconD,sourceFam)

    def addLocusFamily(self, lfO):
        '''Add a LocusFamily. Assumes initializeFamily has already been called
to create the corresponding family.
        '''
        self.locusFamiliesD[lfO.locusFamNum] =  lfO
        self.familiesD[lfO.famNum].addLocusFamily(lfO)

    def getLocusFamily(self,locusFamNum):
        return self.locusFamiliesD[locusFamNum]

    def getFamily(self,famNum):
        return self.familiesD[famNum]

    def iterLocusFamilies(self):
        '''Iterate over all LocusFamilies in order of locusFamNum.
        '''
        maxLocusFamNum = max(self.locusFamiliesD.keys())
        for i in range(maxLocusFamNum+1):
            if i in self.locusFamiliesD:
                yield self.locusFamiliesD[i]
        
    def iterFamilies(self):
        '''Iterate over all Families in order of famNum.'''
        maxFamNum = max(self.familiesD.keys())
        for i in range(maxFamNum+1):
            if i in self.familiesD:
                yield self.familiesD[i]
                
    def getAllGenes(self):
        '''Collect all the genes present in the LocusFamily objects belonging
to this instance of the Families class. Return as a set.
        '''
        allGenesS=set()
        for lfO in self.iterLocusFamilies():
            allGenesS.update(set(lfO.iterGenes()))
        return allGenesS

    def getNumFamilies(self):
        return len(self.familiesD)
    
    def getNumLocusFamilies(self):
        return len(self.locusFamiliesD)
    
    def __repr__(self):
        return "<Families object--"+str(len(self.familiesD))+" Families, "+str(len(self.locusFamiliesD))+" LocusFamilies>"


# temp

def convertReconBranchToNode(brReconD):
    '''Take a reconciliation dictionary in the branch based format. (The
format produced by the dtlor dp function, which represents the
placement of gene tree branches on the species tree). Return a
dictionary where reconciliation events are explicitly characterized as
falling at nodes or branches in the gene and species trees. This output dictionary contains two subdictionaries. One is event based (outD['event']), where keys are of the form

(geneNB,'n/b',speciesNB,'n/b',eventType)

The values include locus numbers on either side of the event and keys for the child events.

Output also contains a second dictionary (outD['species']) where keys are 2-tuples of speciesNB and 'n/b'. And values are (eventKey,locusAtTop,locusAtBottom,cm1Key,cm2Key).
    '''
    rootKeyInputFormat = ''
    for key in brReconD:
        if key[0] == 'root':
            rootKeyInputFormat = key
    
    eventsL = getEventsL(brReconD,rootKeyInputFormat)

    # make two dicts
    rootKey = eventsL[0][0]
    eventD = {}
    speciesD = {}
    for eventKey,locusAtTop,locusAtBottom,cm1Key,cm2Key in eventsL:
        eventD[eventKey] = locusAtTop,locusAtBottom,cm1Key,cm2Key

        speciesNBkey = eventKey[2:4]
        if speciesNBkey in speciesD:
            speciesD[speciesNBkey].append((eventKey,locusAtTop,locusAtBottom,cm1Key,cm2Key))
        else:
            speciesD[speciesNBkey] = [(eventKey,locusAtTop,locusAtBottom,cm1Key,cm2Key)]

    outD = {}
    outD['event'] = eventD
    outD['species'] = speciesD

    return outD
            

def getEventsL(brReconD,key):
    '''Traverse the reconciliation dict in preorder starting at key.'''

    if key[0] == None:
        return [] # reached tip of gene tree
    else:
        _, child1Mapping, child2Mapping = brReconD[key]
        eventT = getOneEvent(brReconD,key)
        lL = getEventsL(brReconD,child1Mapping) # left
        rL = getEventsL(brReconD,child2Mapping) # right
        return [eventT] + lL + rL
        
def getOneEvent(brReconD,key):
    '''Given a key to the branch based reconciliation dictionary, return
all the info in key and value as separate items.'''

    geneNB,speciesNB,locusAtTop,locusAtBottom = key
    eventType, child1Mapping, child2Mapping = brReconD[key]

    outKey = determineNodeOrBranch(geneNB,speciesNB,eventType)
    
    # get event types of children
    cm1GeneNB,cm1spNB,cm1LocAtTop,cm1LocAtBot = child1Mapping # cm1 is childMapping1
    if cm1GeneNB == None:
        cm1Event = None
    else:
        cm1Event = brReconD[child1Mapping][0]

    cm1Key = determineNodeOrBranch(cm1GeneNB,cm1spNB,cm1Event)
        
    cm2GeneNB,cm2spNB,cm2LocAtTop,cm2LocAtBot = child2Mapping # cm2 is childMapping2
    if cm2GeneNB == None:
        cm2Event = None
    else:
        cm2Event = brReconD[child2Mapping][0]

    cm2Key = determineNodeOrBranch(cm2GeneNB,cm2spNB,cm2Event)

    eventT = (outKey,locusAtTop,locusAtBottom,cm1Key,cm2Key)
    
    return eventT

def determineNodeOrBranch(geneNB,speciesNB,eventType):
    '''Determine whether the gene and species locations are on nodes or
branches based on the event type.'''
    
    if eventType == 'D':
        # implies geneNB is a node, speciesNB is a branch
        return (geneNB,'n',speciesNB,'b',eventType)
    elif eventType == 'T':
        # implies geneNB is a node, speciesNB is a branch
        return (geneNB,'n',speciesNB,'b',eventType)
    elif eventType == 'L':
        # implies geneNB is a branch, speciesNB is a node
        return (geneNB,'b',speciesNB,'n',eventType)
    elif eventType == 'O':
        # implies geneNB is a branch, speciesNB is a branch
        return (geneNB,'b',speciesNB,'b',eventType)
    elif eventType == 'R':
        # implies geneNB is a branch, speciesNB is a branch
        return (geneNB,'b',speciesNB,'b',eventType)
    elif eventType == 'C':
        # implies geneNB is a node, speciesNB is a node
        return (geneNB,'n',speciesNB,'n',eventType)
    elif eventType == 'S':
        # implies geneNB is a node, speciesNB is a node
        return (geneNB,'n',speciesNB,'n',eventType)
        # what is S?
    else:
        return # should never get here
