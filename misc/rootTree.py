import sys
from Bio import Phylo

def rootTree(tree,outGroupTaxaL):
    '''Given a bio python tree, root it. outGroupTaxaL should consist of
all the species in the desired outgroup clade.'''

    tree.root_with_outgroup(*outGroupTaxaL)
    return tree 

        
if __name__ == "__main__":
    
    utreeFN = sys.argv[1]
    rtreeFN = sys.argv[2]
    outGroupTaxaL = sys.argv[3:]
    
    utree = Phylo.read(utreeFN, 'newick')
    
    rtree = rootTree(utree,outGroupTaxaL)
    
    Phylo.write(rtree, rtreeFN, 'newick')
