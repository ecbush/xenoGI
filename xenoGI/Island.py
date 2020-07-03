
class LocusIsland:

    def __init__(self, idnum, mrca, locusFamilyL):
        '''Create a LocusIsland object.'''
        self.id = idnum
        self.mrca = mrca
        self.locusFamilyL=locusFamilyL # list of locus family numbers in the island, in order
            
    def __repr__(self):
        return "<id:"+str(self.id)+", mrca:"+str(self.mrca)+", locusFamilyL:"+str(self.locusFamilyL)+">"

    def __len__(self):
        return len(self.locusFamilyL)
    
    def fileStr(self):
        '''Return a string which can be used for saving a island compactly in a
file.'''
        return str(self.id)+"\t"+self.mrca+"\t"+",".join(map(str,self.locusFamilyL))
        
    def merge(self,other,orientation):
        '''Merge island other into self. The argument orientation tells us
which orientation to combine the families of self and other in. The
meaning of the different values is determined by the cases in the
score function.
        '''
        if orientation == 0:
            newFamilyL= self.locusFamilyL + other.locusFamilyL
        elif orientation == 1:
            newFamilyL= self.locusFamilyL + other.locusFamilyL[::-1]
        elif orientation == 2:
            newFamilyL= other.locusFamilyL[::-1] + self.locusFamilyL
        elif orientation == 3:
            newFamilyL= other.locusFamilyL + self.locusFamilyL
        self.locusFamilyL=newFamilyL

    def iterLocusFamilies(self,familiesO):
        '''Iterates, yielding LocusFamily objects.'''
        for locusFamNum in self.locusFamilyL:
            yield familiesO.getLocusFamily(locusFamNum)

    def iterGenes(self,familiesO):
        '''Iterate though all the genes in this LocusIsland.'''
        for lfO in self.iterLocusFamilies(familiesO):
            for gene in lfO.iterGenes():
                yield gene

    def getLocusFamilyOriginStr(self,familiesO,rootFocalClade):
        '''Return a string with the origin for each locus family in the island
(C,X,R).'''
        L = []
        for lfO in self.iterLocusFamilies(familiesO):
            L.append(lfO.origin(familiesO,rootFocalClade))
        return "".join(L)
        
def str2Island(islandStr):
    '''Given a island string (e.g. produced by the fileStr method) parse to produce island.'''
    L=islandStr.split('\t')
    id=int(L[0])
    mrca=L[1]
    locusFamilyL=[int(x) for x in L[2].split(',')]
    return LocusIsland(id,mrca,locusFamilyL)
