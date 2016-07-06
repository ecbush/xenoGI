
class Group:

    def __init__(self, idnum, familyL, mrca=None, familyStrainT=None, tree=None):
        '''New objects must always be provided with an idNum and a familyL. They can either be provided directly with an mrca, or with a familyStrainT and a tree.'''
        self.id = idnum
        self.familyL=familyL # list of family numbers in the group, in order

        if not (mrca!=None or (familyStrainT!=None and tree!=None)):
            raise ValueError('Wrong inputs provided for Group instance. Must provide either mrca, or both familyStrainT and tree.')
            
        if mrca==None:
            # familyL[0] is the identifier for a single family,
            # familyStrainT is the data on gene content of all families,
            # and tree.
            genesInFamilyT=familyStrainT[familyL[0]]
            self.mrca = self.getMRCA(genesInFamilyT,tree)
        else:
            self.mrca = mrca

            
    def __repr__(self):
        return "<id:"+str(self.id)+", mrca:"+str(self.mrca)+", familyL:"+str(self.familyL)+">"

    def fileStr(self):
        '''Return a string which can be used for saving a group compactly in a
file.'''
        return str(self.id)+", "+str(self.mrca)+", "+" ".join(map(str,self.familyL))
        
    def getMRCA(self,genesInFamilyT,tree):
        '''Returns the most recent common ancestor node number for all genes
in the family.
        '''
        if tree[1]==():
            if genesInFamilyT[tree[0]][0]==0: return -1
            else: return tree[0]
        else:
            left=self.getMRCA(genesInFamilyT,tree[1])
            right=self.getMRCA(genesInFamilyT,tree[2])
            if left>-1 and right>-1:
                return tree[0]
            elif left>-1:
                return left
            else: return right

    def merge(self,other,orientation):
        '''Merge group other into self. The argument orientation tells us
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
            
def str2Group(groupStr):
    '''Given a group string (e.g. produced by the fileStr method) parse to produce a group.'''
    L=groupStr.split(', ')
    id=int(L[0])
    mrca=int(L[1])
    familyL=[int(x) for x in L[2].split()]
    return Group(id,familyL,mrca=mrca)

# need to adjust init. it can take mrca is input. or tree and familyStrainT.
