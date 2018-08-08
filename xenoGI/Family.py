
class LocusFamily:
    def __init__(self, famNum, locusFamNum, lfMrca, genesL):
        '''Initialize a LocusFamily object with a family number, a LocusFamily
number, and the mrca for this locus family (which may differ from the
mrca for its family.'''
        self.famNum = famNum
        self.locusFamNum = locusFamNum
        self.lfMrca = lfMrca
        self.genesL = genesL

    def addGene(self, gene):
        '''Add a gene to our set of genes. This should be the numerical
representation of a gene.'''
        self.genesL.append(gene)

    def getGeneNums(self):
        return self.genesL

    def getStr(self,strainNum2StrD,geneNames,sep):
        '''Return a string representation of a single LocusFamily. Separator
between elements given by sep. Elements are: locusFamNum lfMrca gene1
gene2...
        '''
        outL=[str(self.locusFamNum),strainNum2StrD[self.lfMrca]]
        
        for geneNum in self.genesL:
            outL.append(geneNames.numToName(geneNum))

        return sep.join(outL)
    
    def fileStr(self,strainNum2StrD,geneNames):
        '''Return a string representation of a single LocusFamily. Format is
        comma separated: 
        '''
        return self.getStr(strainNum2StrD,geneNames,',')

    def __repr__(self):
        '''String representation of a LocusFamily, for display purposes.'''
        return "<lf "+str(self.locusFamNum)+">"
        
        
class Family:
    def __init__(self,famNum,mrca,seedPairL=None):
        '''Initialize an object of class Family.'''

        self.famNum = famNum
        self.mrca = mrca
        self.seedPairL = seedPairL # original seeds from PHiGS. This will be
                                   # None if single gene family
        self.locusFamiliesL = []   # will contain locusFamily objects for this family

    def addLocusFamily(self,lfO):
        self.locusFamiliesL.append(lfO)

    def getLocusFamilies(self):
        return self.locusFamiliesL

    def getGeneNums(self):
        '''Get all the genes associated with this family in numerical form.'''
        famGenesL=[]
        for lfO in self.getLocusFamilies():
            for gene in lfO.getGeneNums():
                famGenesL.append(gene)
        return famGenesL

    def getOutsideConnections(self,scoresO):
        '''Given a score object, return a set of all outside genes with
connections to this family.'''
        allGenesInFamL = self.getGeneNums()
        otherGenesS=set()
        for geneNum in allGenesInFamL:
            for otherGene in scoresO.getConnectionsGene(geneNum):
                if not otherGene in allGenesInFamL:
                    otherGenesS.add(otherGene)
        return otherGenesS
    
    def fileStr(self,strainNum2StrD,geneNames):
        '''Return string representation of single family. Format is: famNum <tab> 
        mrca <tab> seedG1 <tab> seedG2 <tab> locusFamNum1,locusFamGenes <tab>
        locusFamNum2,locusFamGenes...
        The LocusFamily object representations are comma separated.
        '''
        
        outL =[str(self.famNum)]
        outL.append(strainNum2StrD[self.mrca])

        if self.seedPairL == None:
            outL.extend(["-","-"])
        else:
            for seed in self.seedPairL:
                outL.append(geneNames.numToName(seed))

        for lfO in self.locusFamiliesL:
            outL.append(lfO.fileStr(strainNum2StrD,geneNames))
                
        return "\t".join(outL)

    def __repr__(self):
        return "<fam:"+str(self.famNum)+">"
    

class Families:

    def __init__(self,tree):
        '''Initialize an object of class Families.'''

        self.tree = tree
        self.locusFamiliesD = {}
        self.familiesD = {} 

        ## locusFamiliesD has key locusFamNum and value a LocusFamily
        ## object. familiesD has key famNum and value
        ## [mrca,seedG1,seedG2,[List of lf numbers]]

        
    def initializeFamily(self,famNum,mrca,seedPairL=None):
        '''Set up an entry for family famNum.'''
        # the seed genes are the original PHiGs seed.

        self.familiesD[famNum] = Family(famNum,mrca,seedPairL=seedPairL)

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
            allGenesS.update(lfO.genesL)
        return allGenesS

    def __repr__(self):
        return "<Families object--"+str(len(self.familiesD))+" Families, "+str(len(self.locusFamiliesD))+" LocusFamilies>"
