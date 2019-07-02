
class LocusFamily:
    def __init__(self, famNum, locusFamNum, lfMrca):
        '''Initialize a LocusFamily object with a family number, a LocusFamily
number, and the mrca for this locus family (which may differ from the
mrca for its family.'''
        self.famNum = famNum
        self.locusFamNum = locusFamNum
        self.lfMrca = lfMrca
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
    
    def getStr(self,genesO,sep):
        '''Return a string representation of a single LocusFamily. Separator
between elements given by sep. Elements are: locusFamNum lfMrca gene1
gene2...
        '''
        outL=[str(self.locusFamNum),self.lfMrca]
        
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
    
    def fileStr(self,genesO):
        '''Return string representation of single family. Format is: famNum <tab> 
        mrca <tab> seedG1 <tab> seedG2 <tab> locusFamNum1,locusFamGenes <tab>
        locusFamNum2,locusFamGenes...
        The LocusFamily object representations are comma separated.
        '''
        
        outL =[str(self.famNum),self.mrca]

        if self.seedPairL == None:
            outL.extend(["-","-"])
        else:
            for seed in self.seedPairL:
                outL.append(genesO.numToName(seed))

        for lfO in self.locusFamiliesL:
            outL.append(lfO.fileStr(genesO))
                
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
        ## a Family object.

        
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
            allGenesS.update(set(lfO.iterGenes()))
        return allGenesS

    def __repr__(self):
        return "<Families object--"+str(len(self.familiesD))+" Families, "+str(len(self.locusFamiliesD))+" LocusFamilies>"
