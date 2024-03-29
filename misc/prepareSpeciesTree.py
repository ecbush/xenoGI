import sys,os
from Bio import Phylo
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import trees
from xenoGI import Tree

if __name__ == "__main__":
    
    inTreeFN = sys.argv[1]
    outTreeFN = sys.argv[2]
    outGroupTaxaL = sys.argv[3:]
    
    bpTree = Phylo.read(inTreeFN, 'newick')
    bpTree = trees.prepareTree(bpTree,outGroupTaxaL)
    
    speciesRtreeO = Tree.Rtree()
    nodeConnectD,branchLenD = speciesRtreeO.__bioPhyloToNodeConnectD__(bpTree)    
    speciesRtreeO.populateAttributes(nodeConnectD,bpTree.root.name)
    
    # write final tree. rooted, with named internal nodes and no
    # branch lengths
    with open(outTreeFN,'w') as f:
        f.write(speciesRtreeO.toNewickStr()+";")
