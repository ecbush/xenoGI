import sys
from Bio import Phylo

# assumes a rooted tree

def stripBranchLen(clade):
    '''Remove branch lengths from clade. Modifies in place.'''
    clade.branch_length = 0
    if clade.is_terminal():
        return
    else:
        stripBranchLen(clade[0])
        stripBranchLen(clade[1])
        return

if __name__ == "__main__":

    brTreeFN = sys.argv[1]
    noBrTreeFN = sys.argv[2]

    tree = Phylo.read(brTreeFN, 'newick')
    
    stripBranchLen(tree.root)
    Phylo.write(tree, noBrTreeFN, 'newick',format_branch_length='%d')

    # need one more hacky step to take in string and get rid of :0.
