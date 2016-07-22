
class Family:

    def __setattr__(self, *args):
         raise TypeError("can't modify immutable instance")
    __delattr__ = __setattr__

    def __init__(self, idnum, mrca, genesL,numNodesInTree,geneName2NumD,geneName2StrainNumD):
        '''Given an family ID number, a most recent common ancestor, and a
list of genes (in string form) create an immutable family object. numNodesInTree
is the number of nodes found in the tree (the full tree) we're working
with for the project.
        '''

        super(Family, self).__setattr__('id', idnum)
        super(Family, self).__setattr__('mrca', mrca)
        
        # create the family gene tuple, This has indexes corresponding
        # to nodes on the tree. At each position we have another tuple
        # (gene count, (tuple of genes))

        familyL = [[0,[]] for j in range(numNodesInTree)]
        for gene in genesL:
            geneNum=geneName2NumD[gene]
            strainNum=geneName2StrainNumD[gene]
            familyL[strainNum][0]+=1
            familyL[strainNum][1].append(geneNum)

        # tuple-ize
        newFamilyL=[]
        for ct,L in familyL:
            newFamilyL.append((ct,tuple(L)))
        

        super(Family, self).__setattr__('famGeneT', tuple(newFamilyL))

    def __repr__(self):
        return "Family: "+str(self.id)+' '+str(self.mrca)+' '+str(self.famGeneT)
