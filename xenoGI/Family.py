import sys
from . import trees,new_DTLOR_DP

class LocusFamily:
    def __init__(self, famNum, locusFamNum, lfMrca,locusNum=None,reconRootKey=None):
        '''Initialize a LocusFamily object with a family number, a LocusFamily
number, and the mrca for this locus family (which may differ from the
mrca for its family.'''
        self.famNum = famNum
        self.locusFamNum = locusFamNum # a unique number for each lf
        self.lfMrca = lfMrca # on species tree
        self.locusNum = locusNum # indicates location, several lf may share
        self.geneD = {}
        # root branch of locus family in family recon. (geneTreeLoc,geneTreeNB)
        self.reconRootKey = reconRootKey

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

    def origin(self,familiesO,rootFocalClade):
        '''Determine the event at the origin of this locus family. Returns "C", "X" or
        or "R". If reconRootKey is None, then no reconcilation has
        been done, and we return " ".

        '''
        if self.reconRootKey == None:
            return " "
        else:
            fam = familiesO.getFamily(self.famNum)
            eventTypeL=[]
            for eventValue in fam.reconD[self.reconRootKey]:
                eventType = eventValue[0]
                if eventType == 'R' or eventType == 'O':
                    eventTypeL.append(eventType)
            assert len(eventTypeL) < 2, "Should be at most one O or one R on branch."

            if eventTypeL[0] == 'R':
                orig = 'R'
            else:
                orig = fam.origin(familiesO.speciesRtreeO,rootFocalClade)
            
            return orig

    def printReconByGeneTree(self,familiesO,fileF=sys.stdout):
        '''Print a text summary of the reconciliation corresponding to this
locus family.'''
        if self.reconRootKey == None:
            print(None)
        else:
            fam = familiesO.getFamily(self.famNum)
            rootOfLfSubtree = self.reconRootKey[0]
            lfSubtree = trees.subtree(fam.geneRtree,rootOfLfSubtree)
            fam.printReconByGeneTreeHelper(lfSubtree,0,fileF)
    
    def getStr(self,genesO,sep):
        '''Return a string representation of a single LocusFamily. Separator
between elements given by sep. Elements are: locusFamNum lfMrca gene1
gene2...
        '''
        reconRootKeyStr = ("_".join(map(str,self.reconRootKey)) if self.reconRootKey!=None else str(self.reconRootKey))
        outL=[str(self.locusFamNum),self.lfMrca,str(self.locusNum),reconRootKeyStr]
        
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
    def __init__(self,famNum,mrca,geneTreeO=None,reconD=None):
        '''Base class to be inherited by initialFamily and originFamily.'''

        self.famNum = famNum
        self.mrca = mrca
        self.locusFamiliesL = []   # will contain locusFamily objects for this family
        self.geneTreeO = geneTreeO # gene tree object
        self.reconD = reconD       # dict with reconciliation graph
        
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

    def addGeneTree(self,geneTreeO):
        '''Given tree object geneTreeO, store as attribute.'''
        self.geneTreeO = geneTreeO

    def addReconciliation(self,rD):
        '''Given a dictionary rD storing reconciliation values, store as
attribute.
        '''
        self.reconD = rD 

    def fileStr(self,genesO):
        '''Return string representation of single family. Format is: famNum
        <tab> geneTreeO <tab> reconciliation <tab> mrca <tab>
        locusFamNum1,locusFamGenes <tab> locusFamNum2,locusFamGenes...
        The LocusFamily object representations are comma separated.
        The geneTreeO is string version of an Rtree object.
        reconciliation is a string of the dict.
        '''
        # family number and mrca
        outL =[str(self.famNum),self.mrca]

        # geneTreeO
        if self.geneTreeO != None:
            outL.append(self.geneTreeO.fileStr())
        else:
            outL.append("None")

        # reconciliation
        if self.reconD != None:
            outL.append(str(self.reconD))
        else:
            outL.append("None")

        # dtlorCost or sourceFam
        if hasattr(self,"dtlorCost"):
            val = self.dtlorCost
        else:
            val = self.sourceFam
        
        if val != None:
            outL.append(str(val))
        else:
            outL.append("None")
            
        # locus families
        for lfO in self.locusFamiliesL:
            outL.append(lfO.fileStr(genesO))

        return "\t".join(outL)
    
class initialFamily(Family):

    def __init__(self,famNum,mrca,geneTreeO=None,reconD=None):
        '''Initialize an object of class initialFamily.'''
        super().__init__(famNum,mrca,geneTreeO,reconD)
        self.dtlorCost = None # cost from dtlor alg
        
        # for initial families:
        # geneTreeO is presently a rooted tree. Perhaps in future unrooted.
        # reconD is the raw DTLOR output, a graph dict.

    def countMPRs(self):
        '''How many most parsimonious reconcilations are there with the same
optimal cost?

        '''
        return new_DTLOR_DP.count_MPRs(self.reconD)[(new_DTLOR_DP.NodeType.ROOT,)]
        
    def getRandomOriginMprReconD(self,preOrderT):
        '''From the reconciliation graph (dtlor output), get a randomly chosen
median mpr, convert to our node based format and return as a dict.

        '''
        ## support funcs

        def getEventTypeStr(eventTypeO):
            '''Convert from the object representation (from dtlor) of an event to
a string.
            '''
            if eventTypeO == new_DTLOR_DP.NodeType.DUPLICATION:
                return "D"
            elif eventTypeO == new_DTLOR_DP.NodeType.TRANSFER:
                return "T"
            elif eventTypeO == new_DTLOR_DP.NodeType.LOSS:
                return "L"
            elif eventTypeO == new_DTLOR_DP.NodeType.ORIGIN:
                return "O"
            elif eventTypeO == new_DTLOR_DP.NodeType.REARRANGEMENT:
                return "R"
            elif eventTypeO == new_DTLOR_DP.NodeType.COSPECIATION:
                return "S"
            else:
                # should never get here
                raise ValueError("Event is not among the D T L O R options.")

        def getSpeciesNodeClosestToTips(speciesTreeLocL,postOrderT):
            for stLoc in postOrderT:
                if stLoc in speciesTreeLocL:
                    return stLoc
            # should never get here
            raise ValueError("Species node not in list of postorder nodes from species tree.")

        def determineDTLNodeOrBranch(eventTypeO):
            '''Determine whether the gene and species locations are on nodes or
        branches based on the event type.'''
            if eventTypeO == new_DTLOR_DP.NodeType.DUPLICATION:
                # implies geneTreeLoc is a node, speciesTreeLoc is a branch
                return ('n','b')
            elif eventTypeO == new_DTLOR_DP.NodeType.TRANSFER:
                # implies geneTreeLoc is a node, speciesTreeLoc is a branch
                return ('n','b')
            elif eventTypeO == new_DTLOR_DP.NodeType.LOSS:
                # implies geneTreeLoc is a branch, speciesTreeLoc is a node
                return ('b','n')
            elif eventTypeO == new_DTLOR_DP.NodeType.COSPECIATION:
                # implies geneTreeLoc is a node, speciesTreeLoc is a node
                return ('n','n')
            else:
                return # should never get here
        
        ## main code
        postOrderT = preOrderT[::-1]
        
        # get event list
        eventG = new_DTLOR_DP.build_event_median_graph(self.reconD)
        mpr = new_DTLOR_DP.find_MPR(eventG) # one randomly chosen median mpr
        eventL = new_DTLOR_DP.get_events(mpr)

        # create dict with keys (geneTreeLoc,'n/b')
        reconD = {}
        for eventT in eventL:
            if eventT[0] in (new_DTLOR_DP.NodeType.ORIGIN,new_DTLOR_DP.NodeType.REARRANGEMENT):
                (eventTypeO, geneTreeLoc, location, speciesTreeLocL) = eventT
                eventKey = (geneTreeLoc,'b')
                eventType = getEventTypeStr(eventTypeO)
                speciesTreeLoc = getSpeciesNodeClosestToTips(speciesTreeLocL,postOrderT)
                eventValue = (eventType,speciesTreeLoc,'b',location)
            elif eventT[0] in (new_DTLOR_DP.NodeType.DUPLICATION,new_DTLOR_DP.NodeType.TRANSFER,new_DTLOR_DP.NodeType.LOSS,new_DTLOR_DP.NodeType.COSPECIATION):
                (eventTypeO, geneTreeLoc, speciesTreeLoc) = eventT
                geneTreeNB,speciesTreeNB = determineDTLNodeOrBranch(eventTypeO)
                eventKey = (geneTreeLoc,geneTreeNB)
                eventType = getEventTypeStr(eventTypeO)
                eventValue = (eventType,speciesTreeLoc,speciesTreeNB,None)
                
            # add it
            if eventKey in reconD:
                reconD[eventKey].append(eventValue)
            else:
                reconD[eventKey] = [eventValue]
                
        return reconD

    def __repr__(self):
        return "<ifam:"+str(self.famNum)+">"

class originFamily(Family):

    def __init__(self,famNum,mrca,geneTreeO=None,reconD=None):
        '''Initialize an object of class originFamily.'''
        super().__init__(famNum,mrca,geneTreeO,reconD)

        # for origin families:
        # geneTreeO is a rooted tree.
        # reconD is a converted version of a single DTLOR MPR

        self.sourceFam = None   # ifam num that originFam came from

        # attributes not saved to file and must be recreated when needed
        self.geneHistoryD = None

    def __repr__(self):
        return "<ofam:"+str(self.famNum)+">"

    def createGeneHistoryD(self):
        '''Create a dict of gene histories. Keyed by tip gene. Value is a
string giving history of gene. We convert O events into either C
(core, happened at or before the rootFocalClade) or X (xeno hgt)
events.
        '''
        # support funcs
        def createGeneHistoryStringL(geneRtreeO,node):
            '''Actual recursive function for making strings of gene
histories. Returns a list of tuples. Each tuple contains
(tipGene,historyStr).

            '''
            # get sequence for this branch and node
            hisStr = ""
            for nbKey in [(node,'b'),(node,'n')]:
                if nbKey in self.reconD:
                    for event,stLoc,stNB,synLocus in self.reconD[nbKey]:
                        if event in 'LM':
                            pass
                        else:
                            hisStr+=event
                            
            if geneRtreeO.isLeaf(node):
                return [(node,hisStr)]
            else:
                tempL = []
                for child in geneRtreeO.children(node):
                    partialTempL = createGeneHistoryStringL(geneRtreeO,child)
                    tempL.extend(partialTempL)
                    
                # loop and add current hisStr to events in these
                outL = []
                for geneNum,hisStrRestOfTree in tempL:
                    outL.append((geneNum,hisStr+hisStrRestOfTree))

                return outL
        # end support funcs

        # get history strings
        hisStrL = createGeneHistoryStringL(self.geneTreeO,self.geneTreeO.rootNode)

        # need to only have string back to rootFocalClade in species tree. cut it off before that. do this by keeping the species Tree placements.
        
        # put in dict
        self.geneHistoryD = {}
        for geneNum,hisStr in hisStrL:
            self.geneHistoryD[geneNum] = hisStr

    def getGeneHistoryStr(self,geneNum):
        '''Return a gene history string. Possible characters are:

        D - duplication
        T - transfer (hgt within the species tree)
        O - origin event (entry into species tree)
        R - rearrangment event
        S - cospeciation event

        Note, loss (L) events won't show up in these. And we don't
        print the cotermination (M) events.
        '''

        if self.geneHistoryD != None:
            # geneHistoryD already made
            return self.geneHistoryD[geneNum]    
        elif self.reconD != None:
            # geneHistoryD not created yet and reconciliation is present
            self.createGeneHistoryD()
            return self.geneHistoryD[geneNum]    
        else:
            # this family has no recon, so we can't caculate a history
            return ""

    def origin(self,speciesRtree,rootFocalClade):
        '''Intepret the O event at the base of this family by determining if
it is a core gene (C) or xeno hgt (X). We assess whether it is core or
not at the base of the rootFocalClade of the species tree. If there is
no reconciliation, return empty string.

        '''
        if self.reconD == None:
            return ""

        # get the origin event out. There should only be one. If more,
        # this is probably an intial family.
        Ocount = 0
        for valL in self.reconD.values():
            for eventType,stLoc,stNB,synLocus in valL:
                if eventType == 'O':
                    OstLoc = stLoc
                    Ocount+=1
        if Ocount != 1:
            raise ValueError("There should be exactly one origin event in reconD.")

        # it is C if originated in rfclade branch or ancestors of it
        rfAncT = speciesRtree.ancestors(rootFocalClade) + (rootFocalClade,)
        if OstLoc in rfAncT:
            # core gene
            orig = 'C'
        else:
            orig = 'X'

        return orig
        
    def printReconByGeneTree(self,fileF=sys.stdout):
        '''Print a text summary of the reconciliation.'''
        self.printReconByGeneTreeHelper(self.geneTreeO,self.geneTreeO.rootNode,0,fileF)

    def printReconByGeneTreeHelper(self,geneRtreeO,node,level,fileF=sys.stdout):
        '''Actual recursive function to print reconciliation. Single letter
codes for events are as follows:

        D - duplication
        T - transfer (hgt withing the species tree)
        L - loss (deletion so that gene is lost on species tree lineage)
        O - origin (either core or xeno hgt).
        R - rearrangment event
        M - cotermination event (gene tree tip maps onto species tree tip)

        '''

        def printOneKey(printPrefix,reconD,nbKey):
            if nbKey in reconD:
                for eventT in reconD[nbKey]:
                    gt,gtNB = nbKey
                    event,stLoc,stNB,locus = eventT
                    outStr = printPrefix+event+" ("+str(gt)+" "+gtNB+") "+"("+str(stLoc)+" "+stNB+") "+str(locus)
                    print(outStr,file=fileF)

        # recurse over tree, printing out for each branch and node
        levelSpace = "   "*level
        if geneRtreeO.isLeaf(node):
            print(levelSpace+"- "+node,file=fileF)
            printOneKey(levelSpace+"  ",self.reconD,(node,'b'))
            printOneKey(levelSpace+"  ",self.reconD,(node,'n'))
            return
        else:
            childL = []
            for child in geneRtreeO.children(node):
                childL.append(child)
            childStr = " (children:"+",".join(childL)+")"
            print(levelSpace+"- "+node+childStr,file=fileF)
            printOneKey(levelSpace+"  ",self.reconD,(node,'b'))
            printOneKey(levelSpace+"  ",self.reconD,(node,'n'))
            for child in childL:
                self.printReconByGeneTreeHelper(geneRtreeO,child,level+1,fileF)
    
class Families:

    def __init__(self,speciesRtreeO):
        '''Initialize an object of class Families.'''

        self.speciesRtreeO = speciesRtreeO
        self.locusFamiliesD = {}
        self.familiesD = {} 

        ## locusFamiliesD has key locusFamNum and value a LocusFamily
        ## object. familiesD has key famNum and value
        ## a Family object.
        
    def initializeFamily(self,famNum,mrca,famType,geneTreeO=None,reconD=None,sourceFam=None):
        '''Set up an entry for family famNum.'''

        if famType == "initial":
            self.familiesD[famNum] = initialFamily(famNum,mrca,geneTreeO,reconD)
        else:
            self.familiesD[famNum] = originFamily(famNum,mrca,geneTreeO,reconD)

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

    def addGeneTreeToFamily(self,famNum,geneTreeO):
        self.familiesD[famNum].addGeneTree(geneTreeO)

    def addReconciliationToFamily(self,famNum,reconD):
        self.familiesD[famNum].addReconciliation(reconD)
        
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
