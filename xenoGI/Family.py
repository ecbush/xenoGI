
class LocusFamily:
    def __init__(self, locusFamNum, famNum, lfMrca):
        ''''''

        self.locusFamNum = locusFamNum
        self.famNum = famNum
        self.lfMrca = lfMrca

        self.genesL = []

    def addGene(self, gene):
        '''Add a gene to our set of genes. This should be the numerical
representation of a gene.'''
        self.genesL.append(gene)

    def fileStr(self,strainNum2StrD,geneNames):
        '''Return a string representation of a single LocusFamily. Format is
        comma separated: locusFamNum,gene1,gene2...
        '''

        outL=[str(self.locusFamNum),strainNum2StrD[self.lfMrca]]
        
        for geneNum in self.genesL:
            outL.append(geneNames.numToName(geneNum))

        return ",".join(outL)

        
class Family:
    def __init__(self,famNum,mrca,seedPairL=None):
        '''Initialize an object of class Family.'''

        self.famNum = famNum
        self.mrca = mrca
        self.seedPairL = seedPairL # original seeds from PHiGS. This will be
                                   # None if single gene family
        self.locusFamiliesL = [] # will contain locusFamily IDs for this family

    def addLocusFamily(self,lfO):
        self.locusFamiliesL.append(lfO)
    
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

        for lfO in self.locusfamiliesL:
            outL.append(lfO.fileStr(geneNames))
                
        return "\t".join(outL)


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

    def iterLocusFamily(self):
        '''Iterate over all LocusFamilies in order of locusFamNum.
        '''
        maxLocusFamNum = max(self.locusFamiliesD.keys())
        for i in range(maxLocusFamNum):
            if i in self.locusFamiliesD:
                yield self.locusFamiliesD[i]
        
    def iterFamily(self):
        '''Iterate over all Families in order of famNum.'''
        maxFamNum = max(self.familiesD.keys())
        for i in range(maxFamNum):
            if i in self.familiesD:
                yield self.familiesD[i]
                
    def getAllGenes(self):
        '''Collect all the genes present in the LocusFamily objects belonging
to this instance of the Families class. Return as a set.
        '''
        allGenesS=set()
        for lfO in self.iterLocusFamily():
            allGenesS.update(lfO.genesL)
        return allGenesS
