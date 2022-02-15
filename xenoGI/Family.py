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
            for eventValue in fam.dtlorMprD[self.reconRootKey]:
                eventType = eventValue[0]
                if eventType == 'R' or eventType == 'O':
                    eventTypeL.append(eventType)
            assert len(eventTypeL) < 2, "Should be at most one O or one R on branch."

            if eventTypeL[0] == 'R':
                orig = 'R'
            else:
                orig = fam.origin(familiesO.speciesRtreeO,rootFocalClade)
            
            return orig

    def countEvents(self,familiesO,eventType):
        '''Count the losses in reconciliation associated with this locus family.'''
        fam = familiesO.getFamily(self.famNum)
        # will take reconRootKey as node in counting, since it makes
        # sense to treat O and R events as having happened at very top
        # of their branch.
        if self.reconRootKey == None:
            return None
        else:
            return fam.countEventsBelowNode(self.reconRootKey[0],eventType)

    def dtlorScore(self,familiesO,paramD):
        '''Get dtlor score from reconciliation associated with this locus family.'''
        fam = familiesO.getFamily(self.famNum)
        # will take reconRootKey as node here, since it makes
        # sense to treat O and R events as having happened at very top
        # of their branch.
        if self.reconRootKey == None:
            return None
        else:
            # we need to include the cost for the originating even
            # here because dtlorScoreBelowNode does not include that
            dtlorScore = 0
            if self.origin(familiesO,paramD['rootFocalClade']) in "CX":
                dtlorScore += paramD['originCost']
            elif self.origin(familiesO,paramD['rootFocalClade']) == "R":
                dtlorScore += paramD['rearrangeCost']

            # get scores below
            dtlorScore += fam.dtlorScoreBelowNode(paramD,self.reconRootKey[0])
            
            return dtlorScore
        
    def printReconByGeneTree(self,familiesO,genesO,fileF=sys.stdout):
        '''Print a text summary of the reconciliation corresponding to this
locus family.'''
        if self.reconRootKey == None:
            print(None)
        else:
            fam = familiesO.getFamily(self.famNum)
            rootOfLfSubtree = self.reconRootKey[0]
            lfSubtreeO = fam.geneTreeO.subtree(rootOfLfSubtree)
            fam.printReconByGeneTreeHelper(fam.dtlorMprD,lfSubtreeO,genesO,self.reconRootKey[0],0,fileF)
        #self,dtlorMprD,geneRtreeO,genesO,node,level,fileF=sys.stdout
        
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
    def __init__(self,famNum,mrca,geneTreeO=None,dtlorMprD=None,sourceFam=None):
        '''Base class to be inherited by initialFamily and originFamily.'''

        self.famNum = famNum
        self.mrca = mrca
        self.locusFamiliesL = []   # will contain locusFamily objects for this family
        self.geneTreeO = geneTreeO # gene tree object
        self.dtlorMprD = dtlorMprD # mpr dict
        self.sourceFam = sourceFam # specifies the corresponding fource fam (blast or ifam)
        
    def addLocusFamily(self,lfO):
        self.locusFamiliesL.append(lfO)

    def getLocusFamilies(self):
        return self.locusFamiliesL

    def iterGenes(self):
        '''Iterate through all the genes associated with this family in numerical form.'''
        for lfO in self.getLocusFamilies():
            for gene in lfO.iterGenes():
                yield gene

    def getAllGenes(self):
        '''Collect all the genes present in the LocusFamily objects belonging
to this family. Return as a set.
        '''
        allGenesS=set()
        for lfO in self.getLocusFamilies():
            allGenesS.update(set(lfO.iterGenes()))
        return allGenesS

    def geneCount(self):
        return len(list(self.iterGenes()))
        
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

    def addMprD(self,mprD):
        '''Given a dictionary mprD with an mpr, store as attribute.
        '''
        self.dtlorMprD = mprD 

    def printReconByGeneTreeHelper(self,dtlorMprD,geneRtreeO,genesO,node,level,fileF=sys.stdout):
        '''Actual recursive function to print reconciliation. Single letter
codes for events are as follows:

        D - duplication
        T - transfer (hgt withing the species tree)
        L - loss (deletion so that gene is lost on species tree lineage)
        O - origin (either core or xeno hgt).
        R - rearrangment event
        '''

        def printOneKey(printPrefix,dtlorMprD,nbKey):
            if nbKey in dtlorMprD:
                for eventT in dtlorMprD[nbKey]:
                    gt,gtNB = nbKey
                    event,stLoc,stNB,locus = eventT
                    outStr = printPrefix+event+" ("+str(gt)+" "+gtNB+") --> "+"("+str(stLoc)+" "+stNB+")"
                    if locus != None:
                        outStr+=" synReg:"+str(locus)
                    print(outStr,file=fileF)

        # recurse over tree, printing out for each branch and node
        levelSpace = "   "*level
        if geneRtreeO.isLeaf(node):
            print(levelSpace+"- "+node+" ["+genesO.numToStrainName(int(node))+"]",file=fileF)
            printOneKey(levelSpace+"  ",dtlorMprD,(node,'b'))
            printOneKey(levelSpace+"  ",dtlorMprD,(node,'n'))
            return
        else:
            childL = []
            for child in geneRtreeO.children(node):
                childL.append(child)
            print(levelSpace+"- "+node,file=fileF)
            printOneKey(levelSpace+"  ",dtlorMprD,(node,'b'))
            printOneKey(levelSpace+"  ",dtlorMprD,(node,'n'))
            for child in childL:
                self.printReconByGeneTreeHelper(dtlorMprD,geneRtreeO,genesO,child,level+1,fileF)

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

        # dtlorCost
        if self.dtlorCost != None:
            outL.append(str(self.dtlorCost))
        else:
            outL.append("None")

        # dtlorGraphD
        if self.dtlorGraphD != None:
            outL.append(str(self.dtlorGraphD))
        else:
            outL.append("None")

        # dtlorMprD
        if self.dtlorMprD != None:
            outL.append(str(self.dtlorMprD))
        else:
            outL.append("None")

        # sourceFam
        if self.sourceFam != None:
            outL.append(str(self.sourceFam))
        else:
            outL.append("None")

        # productFamT
        if self.productFamT != None:
            outL.append(str(self.productFamT))
        else:
            outL.append("None")
            
        # locus families
        for lfO in self.locusFamiliesL:
            outL.append(lfO.fileStr(genesO))

        return "\t".join(outL)

class initialFamily(Family):

    def __init__(self,famNum,mrca,geneTreeO=None,dtlorCost=None,dtlorGraphD=None,dtlorMprD=None,sourceFam=None,productFamT=None):
        '''Initialize an object of class initialFamily.'''
        super().__init__(famNum,mrca,geneTreeO,dtlorMprD,sourceFam)
        #  dtlorCost is cost from dtlor alg. A positive integer. In
        #  the case where a family was re-run with permissive origin
        #  costs, we make this value negative.
        self.dtlorCost = dtlorCost
        self.dtlorGraphD = dtlorGraphD
        self.productFamT = productFamT # ids of origin families made from this ifam
        
        # for initial families:
        # geneTreeO is presently a rooted tree. Perhaps in future unrooted.
        # dtlorGraphD is the raw DTLOR output, a graph dict.
        # dtlorMprD is the particular MPR subsequently used
        # sourceFam will be the blast family we came from

    def addGraphD(self,graphD):
        '''Given a dictionary graphD with a dtlor graph, store as attribute.
        '''
        self.dtlorGraphD = graphD 
        
    def productOfams(self):
        '''Return a tuple of the origin families that were made from this
initial family. If that hasn't been done yet, returns None.'''
        return self.productFamT
        
    def countMPRs(self):
        '''How many most parsimonious reconcilations are there with the same
optimal cost?

        '''
        return new_DTLOR_DP.count_MPRs(self.dtlorGraphD)[(new_DTLOR_DP.NodeType.ROOT,)]

    def getMprReconDFromGraph(self,speciesPreOrderT,paramD,isMedian,rand):
        '''From the reconciliation graph (dtlor output), get an mpr, convert
to our node based format. Then return the mpr in original dtlor
format, and also our node based format (both dicts). If isMedian is
True, select a median MPR. If rand is True, this will be randomly
chosen. If rand is False, this will be arbitrarily chosen (and the
same each time if done repeatedly).

        '''

        # get event list
        if isMedian:
            eventG = new_DTLOR_DP.build_event_median_graph(self.dtlorGraphD)
        else:
            eventG = new_DTLOR_DP.build_event_graph(self.dtlorGraphD)

        mprOrigFormatD = new_DTLOR_DP.find_MPR(eventG,rand) # one mpr
        mprNodeFormatD = self.__getMprReconDHelper__(mprOrigFormatD,speciesPreOrderT,paramD)
        
        return mprOrigFormatD,mprNodeFormatD

    def iterMprReconDFromGraph(self,speciesPreOrderT,paramD,isMedian):
        '''Same as getMprReconDFromGraph, but iterates over all MPRs.'''

        # get event list
        if isMedian:
            eventG = new_DTLOR_DP.build_event_median_graph(self.dtlorGraphD)
        else:
            eventG = new_DTLOR_DP.build_event_graph(self.dtlorGraphD)

        for mprOrigFormatD in new_DTLOR_DP.iter_MPRs(eventG):
            # iterate over all
            mprNodeFormatD = self.__getMprReconDHelper__(mprOrigFormatD,speciesPreOrderT,paramD)
        
            yield mprOrigFormatD,mprNodeFormatD
        
    def getMprReconDFromMpr(self,speciesPreOrderT,paramD):
        '''Converts the MPR in self.dtlorMprD to node based format and
returns.
        '''
        return self.__getMprReconDHelper__(self.dtlorMprD,speciesPreOrderT,paramD)

    def printAllPossibleReconsFromGraph(self,speciesPreOrderT,paramD,genesO,fileF=sys.stdout):
        '''Iterate over all possible MPRs in the dtlorGraphD and print.'''

        for mprOrigFormatD,mprNodeFormatD in self.iterMprReconDFromGraph(speciesPreOrderT,paramD,False):
            # do all, not just median mprs
            self.printReconByGeneTreeHelper(mprNodeFormatD,self.geneTreeO,genesO,self.geneTreeO.rootNode,0,fileF)
            print("-------")
            print()
            
    def __getMprReconDHelper__(self,mprOrigFormatD,speciesPreOrderT,paramD):
        '''Takes an MPR in the original format of dtlor and converts to our
node based format.
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
        speciesPostOrderT = speciesPreOrderT[::-1]
                    
        eventL = new_DTLOR_DP.get_events(mprOrigFormatD)

        # create dict with keys (geneTreeLoc,'n/b')
        mprNodeFormatD = {}
        for eventT in eventL:
            if eventT[0] in (new_DTLOR_DP.NodeType.ORIGIN,new_DTLOR_DP.NodeType.REARRANGEMENT):
                (eventTypeO, geneTreeLoc, location, speciesTreeLocL) = eventT
                eventKey = (geneTreeLoc,'b')
                eventType = getEventTypeStr(eventTypeO)
                speciesTreeLoc = getSpeciesNodeClosestToTips(speciesTreeLocL,speciesPostOrderT)
                eventValue = (eventType,speciesTreeLoc,'b',location)
            elif eventT[0] in (new_DTLOR_DP.NodeType.DUPLICATION,new_DTLOR_DP.NodeType.TRANSFER,new_DTLOR_DP.NodeType.LOSS,new_DTLOR_DP.NodeType.COSPECIATION):
                (eventTypeO, geneTreeLoc, speciesTreeLoc) = eventT
                geneTreeNB,speciesTreeNB = determineDTLNodeOrBranch(eventTypeO)
                eventKey = (geneTreeLoc,geneTreeNB)
                eventType = getEventTypeStr(eventTypeO)
                eventValue = (eventType,speciesTreeLoc,speciesTreeNB,None)
                
            # add it
            if eventKey in mprNodeFormatD:
                mprNodeFormatD[eventKey].append(eventValue)
            else:
                mprNodeFormatD[eventKey] = [eventValue]

        # make sure the costs implied by mprNodeFormatD correspond to what
        # dtlor alg gave
        self.__costCheck__(mprNodeFormatD,paramD)
        
        return mprNodeFormatD

    def __costCheck__(self,mprNodeFormatD,paramD):
        '''Sanity check to ensure the events in mprNodeFormatD sum to dtlorCost.'''

        # get costs
        if self.dtlorCost > 0:
            # use normal costs
            D=int(paramD["duplicationCost"])
            T=int(paramD["transferCost"])
            L=int(paramD["lossCost"])
            O=int(paramD["originCost"])
            R=int(paramD["rearrangeCost"])
            dtlorCost = self.dtlorCost
            
        else:
            # neg dtlorCost means we used permissive costs
            D=int(paramD["DTLRcostPermissiveOrigin"])
            T=int(paramD["DTLRcostPermissiveOrigin"])
            L=int(paramD["DTLRcostPermissiveOrigin"])
            O=int(paramD["originCostPermissiveOrigin"])
            R=int(paramD["DTLRcostPermissiveOrigin"])
            dtlorCost = -self.dtlorCost # remove the negative
            
        sm = self.__costSum__(mprNodeFormatD,D,T,L,O,R)
        assert round(sm,3) == round(dtlorCost,3), "Events in this reconcilation don't sum to minCost."

    def __costSum__(self,mprNodeFormatD,D,T,L,O,R):
        '''Sum the costs implied by the recon in mprNodeFormatD.'''
        sm = 0
        for valL in mprNodeFormatD.values():
            for eventType,_,_,_ in valL:
                if eventType == 'D':
                    sm+=D
                elif eventType == 'T':
                    sm+=T
                elif eventType == 'L':
                    sm+=L
                elif eventType == 'O':
                    sm+=O
                elif eventType == 'R':
                    sm+=R
                # S events have no cost
        return sm
    
    def __repr__(self):
        return "<ifam:"+str(self.famNum)+">"

class originFamily(Family):

    def __init__(self,famNum,mrca,geneTreeO=None,dtlorMprD=None,sourceFam=None):
        '''Initialize an object of class originFamily.'''
        super().__init__(famNum,mrca,geneTreeO,dtlorMprD,sourceFam)

        # for origin families:
        # geneTreeO is a rooted tree.
        # dtlorMprD is a converted version of a single DTLOR MPR
        # sourceFam will hold the source ifam for an originFamily

        self.sourceFam = sourceFam # ifam num that originFam came from

        # set to None
        self.dtlorCost = None
        self.dtlorGraphD = None
        self.productFamT = None
        
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
                if nbKey in self.dtlorMprD:
                    for event,stLoc,stNB,synLocus in self.dtlorMprD[nbKey]:
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
        elif self.dtlorMprD != None:
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
        if self.dtlorMprD == None:
            return ""

        # get the origin event out. There should only be one. If more,
        # this is probably an intial family.
        Ocount = 0
        for valL in self.dtlorMprD.values():
            for eventType,stLoc,stNB,synLocus in valL:
                if eventType == 'O':
                    OstLoc = stLoc
                    Ocount+=1
        if Ocount != 1:
            raise ValueError("There should be exactly one origin event in dtlorMprD.")

        # it is C if originated in rfclade branch or ancestors of it
        rfAncT = speciesRtree.ancestors(rootFocalClade)
        if OstLoc in rfAncT:
            # core gene
            orig = 'C'
        else:
            orig = 'X'

        return orig
        
    def printReconByGeneTree(self,genesO,fileF=sys.stdout):
        '''Print a text summary of the reconciliation.'''
        self.printReconByGeneTreeHelper(self.dtlorMprD,self.geneTreeO,genesO,self.geneTreeO.rootNode,0,fileF)

    def countEventsBelowNode(self,geneTreeNode,eventType):
        '''Count how many events of eventType occur in the reconcilation below
geneTreeNode. eventType is a string, D T L O R.'''

        if self.geneTreeO.isLeaf(geneTreeNode):
            return 0
        else:
            eventCtr = 0
            for child in self.geneTreeO.children(geneTreeNode):
                # count on child branch, then recurse
                eventValueL = []
                if (child,'b') in self.dtlorMprD:
                    eventValueL.extend(self.dtlorMprD[(child,'b')])
                if (child,'n') in self.dtlorMprD:
                    eventValueL.extend(self.dtlorMprD[(child,'n')])
                for eventValue in eventValueL:
                    if eventValue[0] == eventType:
                        eventCtr += 1
                eventCtr += self.countEventsBelowNode(child,eventType)
            return eventCtr

    def dtlorScoreBelowNode(self,paramD,geneTreeNode):
        '''Compile dtlor score in the reconcilation below geneTreeNode.
        '''
        if self.geneTreeO.isLeaf(geneTreeNode):
            return 0
        else:
            dtlorScore = 0
            for child in self.geneTreeO.children(geneTreeNode):
                # count on child branch, then recurse
                if (child,'b') in self.dtlorMprD:
                    for eventValue in self.dtlorMprD[(child,'b')]:
                        if eventValue[0] == 'D':
                            dtlorScore += paramD['duplicationCost']
                        elif  eventValue[0] == 'T':
                            dtlorScore += paramD['transferCost']
                        elif  eventValue[0] == 'L':
                            dtlorScore += paramD['lossCost']
                        elif  eventValue[0] == 'O':
                            dtlorScore += paramD['originCost']
                        elif  eventValue[0] == 'R':
                            dtlorScore += paramD['rearrangeCost']
                dtlorScore += self.dtlorScoreBelowNode(paramD,child)
            return dtlorScore
        
    def getNewickGeneTreeWithReconLabels(self,genesO,includeBrLength=False):
        '''Print out a Newick string where the nodes are labeled with
events. Format for annotations are [branch events | node events ]. For
tips on the gene tree, the "node event" is the species it ends up
in.
        '''
        
        # construct nodeLabelD
        nodeLabelD = {}
        for node in self.geneTreeO.preorder():
            nodeStr=""
            # branch events
            if (node,'b') in self.dtlorMprD:
                for eventT in self.dtlorMprD[(node,'b')]:
                    nodeStr += eventT[0]

            # node events
            nodeStr+="|"
            if (node,'n') in self.dtlorMprD:
                for eventT in self.dtlorMprD[(node,'n')]:
                    nodeStr += eventT[0]
                    
            if self.geneTreeO.isLeaf(node):
                nodeStr += genesO.numToStrainName(int(node))
                    
            nodeLabelD[node] = nodeStr

        # put in newick string
        return self.geneTreeO.toNewickStr(includeBrLength=includeBrLength,nodeLabelD=nodeLabelD)

class Families:

    def __init__(self,speciesRtreeO):
        '''Initialize an object of class Families.'''

        self.speciesRtreeO = speciesRtreeO
        self.locusFamiliesD = {}
        self.familiesD = {} 

        ## locusFamiliesD has key locusFamNum and value a LocusFamily
        ## object. familiesD has key famNum and value
        ## a Family object.
        
    def initializeFamily(self,famNum,mrca,famType,geneTreeO=None,dtlorCost=None,dtlorGraphD=None,dtlorMprD=None,sourceFam=None,productFamT=None):
        '''Set up an entry for family famNum.'''

        if famType == "initial":
            self.familiesD[famNum] = initialFamily(famNum,mrca,geneTreeO=geneTreeO,dtlorCost=dtlorCost,dtlorGraphD=dtlorGraphD,dtlorMprD=dtlorMprD,sourceFam=sourceFam,productFamT=productFamT)
        else:
            self.familiesD[famNum] = originFamily(famNum,mrca,geneTreeO=geneTreeO,dtlorMprD=dtlorMprD,sourceFam=sourceFam)

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

    def addMprToFamily(self,famNum,dtlorMprD):
        self.familiesD[famNum].addMprD(dtlorMprD)

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

    def labelHardCore(self,aabrhHardCoreL,typeToLabel):
        '''Run through all families or locusFamiles (specified by typeToLabel)
and create aabrhHardCore attributes in each, indicating which if any
correspond to aabrh hard core families. This attribute can be [], or a
list containing aabrh family numbers.'''
        # funcs
        def matchOneAabrh(aabrhInd,aabrhS,flfItr):
            '''Match single aabrh to objects in flfItr. If match, set atribute of
               that object to True. Return True if match, False if none.
            '''
            for flfO in flfItr():
                S = set(flfO.iterGenes())
                I = S.intersection(aabrhS)
                if len(I) == len(aabrhS):
                    flfO.aabrhHardCore.append(aabrhInd)
                    return True
            return False
        
        # labelHardCore main
        if typeToLabel == 'locusFamily':
            flfItr = self.iterLocusFamilies
        elif typeToLabel == 'Family':
            flfItr = self.iterFamilies

        # create aabrhHardCore attribute, and set to None
        for flfO in flfItr():
            flfO.aabrhHardCore = []

        # now match
        for i,aabrhT in enumerate(aabrhHardCoreL):
            aabrhS = set(aabrhT)
            matchOneAabrh(i,aabrhS,flfItr) # records in object attribute aabrhHardCoreB
            #if not matchOneAabrh(i,aabrhS,flfItr):
            #    print("aabrhS",i,"has no match")
                    
    def __repr__(self):
        return "<Families object--"+str(len(self.familiesD))+" Families, "+str(len(self.locusFamiliesD))+" LocusFamilies>"
