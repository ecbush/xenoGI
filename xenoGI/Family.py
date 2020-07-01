import sys
from . import trees

class LocusFamily:
    def __init__(self, famNum, locusFamNum, lfMrca,locusNum=None,reconRootKey=None):
        '''Initialize a LocusFamily object with a family number, a LocusFamily
number, and the mrca for this locus family (which may differ from the
mrca for its family.'''
        self.famNum = famNum
        self.locusFamNum = locusFamNum # a unique number for each lf
        self.lfMrca = lfMrca # species tree
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

    def origin(self,familiesO):
        '''Determine the event at the origin of this locus family. Returns "O"
        or "R". If reconRootKey is None, then no reconcilation has
        been done, and we return None.

        '''
        if self.reconRootKey == None:
            return None
        else:
            fam = familiesO.getFamily(self.famNum)
            eventTypeL=[]
            for eventValue in fam.reconD[self.reconRootKey]:
                eventType = eventValue[0]
                if eventType == 'R' or eventType == 'O':
                    eventTypeL.append(eventType)
            assert len(eventTypeL) < 2, "Should be at most one O or one R on branch."
            return eventTypeL[0]

    def printReconByGeneTree(self,familiesO,fileF=sys.stdout):
        '''Print a text summary of the reconciliation corresponding to this
locus family.'''
        if self.reconRootKey == None:
            print(None)
        else:
            fam = familiesO.getFamily(self.famNum)
            rootOfLfSubtree = self.reconRootKey[0]
            lfSubtree = trees.subtree(fam.geneTree,rootOfLfSubtree)
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
    def __init__(self,famNum,mrca,geneTree=None,reconD=None,sourceFam=None):
        '''Initialize an object of class Family.'''

        self.famNum = famNum
        self.mrca = mrca
        self.locusFamiliesL = []   # will contain locusFamily objects for this family
        self.geneTree = geneTree   # rooted gene tree object
        self.reconD = reconD # mapping of gene tree onto species tree.
        self.sourceFam = sourceFam # only for originFams. gives ifam num that originFam came from

        # attributes not saved to file and must be recreated when needed
        self.geneHistoryD = None
        
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
        
    def createGeneHistoryD(self,speciesTree,rootFocalClade):
        '''Create a dict of gene histories. Keyed by tip gene. Value is a
string giving history of gene. We convert O events into either C
(core, happened at or before the rootFocalClade) or X (xeno hgt)
events.
        '''
        # support func
        def createGeneHistoryStringL(rfCladeAndAncestorsL,geneTree):
            '''Actual recursive function for making strings of gene
histories. Returns a list of tuples. Each tuple contains
(tipGene,historyStr).

            '''
            # get sequence for this branch and node
            hisStr=''
            for nbKey in [(geneTree[0],'b'),(geneTree[0],'n')]:
                if nbKey in self.reconD:
                    for event,stLoc,stNB,locus in self.reconD[nbKey]:
                        if event == 'L':
                            pass
                        elif event == 'O':
                            if stLoc in rfCladeAndAncestorsL:
                                hisStr+='C' # core gene
                            else:
                                hisStr+='X' # xeno hgt event
                        else:
                            hisStr+=event

            if geneTree[1] == ():
                return [(geneTree[0],hisStr)]
            else:
                lL = createGeneHistoryStringL(rfCladeAndAncestorsL,geneTree[1])
                rL = createGeneHistoryStringL(rfCladeAndAncestorsL,geneTree[2])

                # loop and add current hisStr to events in these
                outL = []
                for geneNum,hisStrRestOfTree in lL+rL:
                    outL.append((geneNum,hisStr+hisStrRestOfTree))

                return outL
        # end support func
            
        # get list of nodes ancestral to and including rootFocalClade
        rfCladeAndAncestorsL = trees.ancestors(speciesTree,rootFocalClade) + [rootFocalClade]

        # get history strings
        hisStrL = createGeneHistoryStringL(rfCladeAndAncestorsL,self.geneTree)

        # need to only have string back to rootFocalClade in species tree. cut it off before that. do this by keeping the speciesTree placements.
        
        # put in dict
        self.geneHistoryD = {}
        for geneNum,hisStr in hisStrL:
            self.geneHistoryD[geneNum] = hisStr

    def getGeneHistoryStr(self,geneNum,speciesTree,rootFocalClade):
        '''Return a gene history string. Possible characters are:

        D - duplication
        T - transfer (hgt within the species tree)
        C - core gene origin event
        X - xeno hgt origin event
        R - rearrangment event

        Note, loss (L) events won't show up in these. And we don't
        print the cotermination (M) events.
        '''

        if self.geneHistoryD != None:
            # geneHistoryD already made
            return self.geneHistoryD[geneNum]    
        elif self.reconD != None:
            # geneHistoryD not created yet and reconciliation is present
            self.createGeneHistoryD(speciesTree,rootFocalClade)
            return self.geneHistoryD[geneNum]    
        else:
            # this family has no recon, so we can't caculate a history
            return ""
        
    def printReconByGeneTree(self,fileF=sys.stdout):
        '''Print a text summary of the reconciliation.'''
        self.printReconByGeneTreeHelper(self.geneTree,0,fileF)

    def printReconByGeneTreeHelper(self,geneTree,level,fileF=sys.stdout):
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
        if geneTree[1] == ():
            print(levelSpace+"- "+str(geneTree[0]),file=fileF)
            printOneKey(levelSpace+"  ",self.reconD,(geneTree[0],'b'))
            printOneKey(levelSpace+"  ",self.reconD,(geneTree[0],'n'))
            return
        else:
            childStr = " (children:"+str(geneTree[1][0])+","+str(geneTree[2][0])+")"
            print(levelSpace+"- "+str(geneTree[0])+childStr,file=fileF)
            printOneKey(levelSpace+"  ",self.reconD,(geneTree[0],'b'))
            printOneKey(levelSpace+"  ",self.reconD,(geneTree[0],'n'))

            self.printReconByGeneTreeHelper(geneTree[1],level+1,fileF)
            self.printReconByGeneTreeHelper(geneTree[2],level+1,fileF)
    
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

    def addGeneTreeToFamily(self,famNum,geneTree):
        self.familiesD[famNum].addGeneTree(geneTree)

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
