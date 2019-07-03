import sys,os
from Bio import Phylo
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import trees


if __name__ == "__main__":
    
    inTreeFN = sys.argv[1]
    outTreeFN = sys.argv[2]
    outGroupTaxaL = sys.argv[3:]
    
    bpTree = Phylo.read(inTreeFN, 'newick')
    
    tupleTree=trees.prepareTree(bpTree,outGroupTaxaL)

    # write final tree. rooted, with named internal nodes and no
    # branch lengths
    trees.writeTreeNoBrLen(tupleTree,outTreeFN)
