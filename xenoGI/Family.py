import sys

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
    def __init__(self,famNum,mrca,geneTree=None,reconD={},sourceFam=None):
        '''Initialize an object of class Family.'''

        self.famNum = famNum
        self.mrca = mrca
        self.locusFamiliesL = []   # will contain locusFamily objects for this family
        self.geneTree = geneTree   # rooted gene tree object
        # dict giving mapping of gene tree onto species tree. May have
        # keys rawBranch, event, species corresponding to different
        # formats. e.g. value of rawBranch is the output format of
        # DTLOR_DP
        self.reconD = reconD
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

    def addReconciliation(self,rD,key):
        '''Given a dictionary rD storing reconciliation values, store as attribute.'''
        self.reconD[key] = rD
        
    def convertReconBranchToNode(self):
        '''Take a reconciliation dictionary in the branch based format. (The
    format produced by the dtlor dp function, which represents the
    placement of gene tree branches on the species tree). Return a
    dictionary where reconciliation events are explicitly
    characterized as falling at nodes or branches in the gene and
    species trees. This output dictionary contains two
    subdictionaries. One is event based (outD['event']), where keys
    are of the form

    (geneNB,'n/b',speciesNB,'n/b',eventType)

    The values (topIsStarB,locusAtBottom,cm1Key,cm2Key) where cm1Key
        is a child mapping.

        '''

        if 'rawBranch' not in self.reconD:
            raise ValueError("Reconciliation dict does not contain reconciliation in raw branch format.")
        
        rootKeyInputFormat = ''
        for key in self.reconD['rawBranch']:
            if key[0] == 'root':
                rootKeyInputFormat = key

        eventsL = self.getEventsL(rootKeyInputFormat)

        # make output dict
        rootKey = eventsL[0][0]
        eventD = {}
        for eventKey,topIsStarB,locusAtBottom,cm1Key,cm2Key in eventsL:
            eventD[eventKey] = topIsStarB,locusAtBottom,cm1Key,cm2Key

        return eventD

    def getEventsL(self,brKey):
        '''Traverse the branch based reconciliation dict in preorder starting
    at brKey.'''

        if brKey[0] == None:
            return [] # reached tip of gene tree
        else:
            _, child1Mapping, child2Mapping = self.reconD['rawBranch'][brKey]
            eventTL = self.parseOneBranchEntry(brKey)
            lL = self.getEventsL(child1Mapping) # left
            rL = self.getEventsL(child2Mapping) # right
            return eventTL + lL + rL

    def parseOneBranchEntry(self,brKey):
        '''Given a key to the branch based reconciliation dictionary, extract
    the event or events (can have multiple when O or R occurs), and return
    as a list of tuples. Each tuples has all the info for one event.

        '''
        eventTL = [] # output
        geneNB,speciesNB,topIsStarB,locusAtBottom = brKey
        eventType, child1Mapping, child2Mapping = self.reconD['rawBranch'][brKey]

        # get the key we'll use for main event (DTLSC)
        eventKey = self.determineNodeOrBranch(geneNB,speciesNB,eventType)


        # temp until the brReconDs have O and R explicitly
        if topIsStarB==True and locusAtBottom!='*':
            # O event
            inducedEventKey = self.determineNodeOrBranch(geneNB,speciesNB,'O')

            eventT = (inducedEventKey,topIsStarB,locusAtBottom,eventKey,None) # points to associated DTS
            eventTL.append(eventT)
        # I'm skipping R, will do when they update the beReconD format.

        # get event types of children
        cm1GeneNB,cm1spNB,cm1LocAtTop,cm1LocAtBot = child1Mapping # cm1 is childMapping1
        if cm1GeneNB == None:
            cm1Event = None
        else:
            cm1Event = self.reconD['rawBranch'][child1Mapping][0]

        cm1Key = self.determineNodeOrBranch(cm1GeneNB,cm1spNB,cm1Event)

        cm2GeneNB,cm2spNB,cm2LocAtTop,cm2LocAtBot = child2Mapping # cm2 is childMapping2
        if cm2GeneNB == None:
            cm2Event = None
        else:
            cm2Event = self.reconD['rawBranch'][child2Mapping][0]

        cm2Key = self.determineNodeOrBranch(cm2GeneNB,cm2spNB,cm2Event)

        eventT = (eventKey,topIsStarB,locusAtBottom,cm1Key,cm2Key)
        eventTL.append(eventT)

        return eventTL

    def determineNodeOrBranch(self,geneNB,speciesNB,eventType):
        '''Determine whether the gene and species locations are on nodes or
    branches based on the event type.'''

        if eventType == 'D': # duplication
            # implies geneNB is a node, speciesNB is a branch
            return (geneNB,'n',speciesNB,'b',eventType)
        elif eventType == 'T': # transfer
            # implies geneNB is a node, speciesNB is a branch
            return (geneNB,'n',speciesNB,'b',eventType)
        elif eventType == 'L': # loss
            # implies geneNB is a branch, speciesNB is a node
            return (geneNB,'b',speciesNB,'n',eventType)
        elif eventType == 'O': # origin
            # implies geneNB is a branch, speciesNB is a branch
            return (geneNB,'b',speciesNB,'b',eventType)
        elif eventType == 'R': # rearrangement
            # implies geneNB is a branch, speciesNB is a branch
            return (geneNB,'b',speciesNB,'b',eventType)
        elif eventType == 'C': # cotermination
            # implies geneNB is a node, speciesNB is a node
            return (geneNB,'n',speciesNB,'n',eventType)
        elif eventType == 'S': # cospeciation
            # implies geneNB is a node, speciesNB is a node
            return (geneNB,'n',speciesNB,'n',eventType)
        else:
            return # should never get here

    def printReconByGeneTree(self,eventD,fileF=sys.stdout):
        ''''''

        # find starting point
        rootkey = ''
        for key in eventD:
            if key[0] == 'root':
                rootkey = key

        for l in self.traverseReconByGeneTree(eventD,rootkey):
            print(l,file=fileF)
            
    def traverseReconByGeneTree(self,eventD,key):
        ''''''

        if key[4] == 'C':
            # cotermination, stop.
            return [(key,eventD[key])]
        else:
            L = [(key,eventD[key])]
            _,_,lkey,rkey = eventD[key]
            lL = self.traverseReconByGeneTree(eventD,lkey)
            
            if key[4] in 'DTS':
                # branching event
                rL = self.traverseReconByGeneTree(eventD,rkey)
            else:
                # is O or R, no event in right slot
                rL = []
            return L + lL + rL
            
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
        
    def initializeFamily(self,famNum,mrca,geneTree=None,reconD={},sourceFam=None):
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
