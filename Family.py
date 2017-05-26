
class Family:

    def __init__(self, idnum, mrca, genesL,numNodesInTree,geneNames,possibleErrorCt):
        '''Given an family ID number, a most recent common ancestor, and a
list of genes create a family object. Handles genes specified either
numerically or by name. numNodesInTree is the number of nodes found in
the tree (the full tree) we're working with for the
project. possibleErrorCt is a count of near miss genes, that either
were added but almost weren't, or were not added but almost were.
        '''

        self.id = idnum
        self.mrca = mrca

        # possibleErrorCt applies to multi-gene families. For single
        # gene families it is None.
        self.possibleErrorCt = possibleErrorCt 
        
        # create the family gene tuple, This has indexes corresponding
        # to nodes on the tree. At each position we have another tuple
        # (gene count, (tuple of genes))

        familyL = [[0,[]] for j in range(numNodesInTree)]
        for gene in genesL:
            # if we're being passed in genes with string names,
            # convert to number
            if type(gene)==int:
                geneNum=gene
            else:
                geneNum=geneNames.nameToNum(gene)

            strainNum=geneNames.numToStrainNum(geneNum)
            familyL[strainNum][0]+=1
            familyL[strainNum][1].append(geneNum)

        # tuple-ize
        newFamilyL=[]
        for ct,L in familyL:
            newFamilyL.append((ct,tuple(L)))
        
        self.famGeneT = tuple(newFamilyL)

    def __repr__(self):
        '''String representation of a family containing family number.'''
        return "<Family: "+str(self.id) + ">"

    def fileStr(self,geneNames,strainNum2StrD):
        '''String representation suitable for writing the family to
file. Genes and mrca are expressed in word form.'''

        outL =[str(self.id)]
        outL.append(strainNum2StrD[self.mrca])
        outL.append(str(self.possibleErrorCt))
        for ct,geneT in self.famGeneT:
            for gene in geneT:
                outL.append(geneNames.numToName(gene))
        return "\t".join(outL)
