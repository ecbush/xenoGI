
class Island:

    def __init__(self, idnum, mrca, familyL):
        '''Create a island object.'''
        self.id = idnum
        self.mrca = mrca
        self.familyL=familyL # list of family numbers in the island, in order
            
    def __repr__(self):
        return "<id:"+str(self.id)+", mrca:"+str(self.mrca)+", familyL:"+str(self.familyL)+">"

    def __len__(self):
        return len(self.familyL)
    
    def fileStr(self,strainNum2StrD):
        '''Return a string which can be used for saving a island compactly in a
file.'''
        return str(self.id)+"\t"+strainNum2StrD[self.mrca]+"\t"+",".join(map(str,self.familyL))
        
    def merge(self,other,orientation):
        '''Merge island other into self. The argument orientation tells us
which orientation to combine the families of self and other in. The
meaning of the different values is determined by the cases in the
score function.
        '''
        if orientation == 0:
            newFamilyL= self.familyL + other.familyL
        elif orientation == 1:
            newFamilyL= self.familyL + other.familyL[::-1]
        elif orientation == 2:
            newFamilyL= other.familyL[::-1] + self.familyL
        elif orientation == 3:
            newFamilyL= other.familyL + self.familyL
        self.familyL=newFamilyL
            
def str2Island(islandStr,strainStr2NumD):
    '''Given a island string (e.g. produced by the fileStr method) parse to produce a island.'''
    L=islandStr.split('\t')
    id=int(L[0])
    mrca=strainStr2NumD[L[1]]
    familyL=[int(x) for x in L[2].split(',')]
    return Island(id,mrca,familyL)
