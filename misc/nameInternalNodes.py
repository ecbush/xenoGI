import sys
from Bio import Phylo

def nameInternalNodes(tree):
    '''Given a bio python tree without named internal nodes, name them.'''
    for j,n in enumerate(tree.get_nonterminals()):
        if n.name != None:
            # this one's already got a name
            raise ValueError("Tree already has at least one named internal node.")
        n.name="i"+str(j)
    return tree 

        
if __name__ == "__main__":
    
    inTreeFN = sys.argv[1]
    outTreeFN = sys.argv[2]
    
    inTree = Phylo.read(inTreeFN, 'newick', rooted=True)
    
    outTree = nameInternalNodes(inTree)
    
    Phylo.write(outTree, outTreeFN, 'newick')
