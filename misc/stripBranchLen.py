import sys,os
from Bio import Phylo
sys.path.append(os.path.join(sys.path[0],'..'))
import trees


def tupleTree2NoBrLenNewick(tree):
    '''Convert a four tuple based tree (root,left,right,branchLen) into a
newick formated string, without branch lengths.'''
    if tree[1]==():
        return str(tree[0])
    else:
        leftStr=tupleTree2NoBrLenNewick(tree[1])
        rightStr=tupleTree2NoBrLenNewick(tree[2])
        return "("+leftStr+","+rightStr+")"+str(tree[0])

def writeTreeNoBrLen(tree,fileName):
    '''Write tree to fileName (in newick format).'''
    f=open(fileName,"w")
    f.write(tupleTree2NoBrLenNewick(tree)+";")
    f.close()

    
def stripBranchLen(tree):
    '''Remove branch lengths from tree. Modifies in place.'''
    for n in tree.get_nonterminals() + tree.get_terminals():
        del n.branch_length
    return tree
        
if __name__ == "__main__":

    inTreeFN = sys.argv[1]
    outTreeFN = sys.argv[2]

    bpTree = Phylo.read(inTreeFN, 'newick', rooted=True)

    trees.checkTree(bpTree)

    #sys.exit()
    
    stringTupleTree=trees.bioPhyloToTupleTree(bpTree)
     
    writeTreeNoBrLen(stringTupleTree,outTreeFN)

