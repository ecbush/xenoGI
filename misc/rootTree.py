import sys
from Bio import Phylo

def rootTree(tree,outGroupTaxaL):
    '''Given a bio python tree, root it. outGroupTaxaL should consist of
all the species in the desired outgroup clade.'''

    tree.root_with_outgroup(*outGroupTaxaL)
    # rooting with outgroups with multiple species seems to leave a
    # trifurcating root.
    if len(tree.root)==2:
        return tree
    else:
        return makeRootBifurcation(tree,outGroupTaxaL)        

def makeRootBifurcation(tree,outGroupTaxaL):
    '''Given a tree that is not bifurcating at root, if it is trifurcating convert to bifurcating. Otherwise thow error.'''

    if len(tree.root)!=3:
        raise ValueError("This tree has a strange number of branches at the root (not 2 or 3).")

    # find the two in outGroupTaxaL, and the one not in it
    inL=[]
    notIn=None
    for cld in tree.root:
        if tipsInOutgroupL(cld,outGroupTaxaL):
            inL.append(cld)
        else:
            notIn = cld

    if len(inL)!=2:
        raise ValueError("There should be two clades of this trifurcating tree which are in outGroupTaxaL, but there are not.")

    # We should take the branch length of notIn and split it up, one
    # half going to notIn, the other half to the new inL clade

    brToSplit=notIn.branch_length
    notIn.branch_length = brToSplit/2
    
    c1=Phylo.Newick.Clade(clades=inL,branch_length=brToSplit/2)
    c2=Phylo.Newick.Clade(clades=[c1,notIn],branch_length=0)
    outTree = Phylo.Newick.Tree(root=c2,rooted=True)

    return outTree
    
def tipsInOutgroupL(clade,outGroupTaxaL):
    '''Return True if all of the tips in clade are in outGroupTaxaL False
otherwise..'''
    for tipName in cladeTipNames(clade):
        if tipName not in outGroupTaxaL:
            return False
    return True
    
    
def cladeTipNames(clade):
    '''Get the names of all the tips in this clade.'''
    L=[]
    for tipC in clade.get_terminals():
        L.append(tipC.name)
    return L
    
if __name__ == "__main__":
    
    utreeFN = sys.argv[1]
    rtreeFN = sys.argv[2]
    outGroupTaxaL = sys.argv[3:]
    
    utree = Phylo.read(utreeFN, 'newick')

    if not utree.is_bifurcating():
        raise ValueError("This tree is not bifurcating. A trifurcation at the root is permitted (to signify it's unrooted) and is not the problem here. There cannot be other trifurcations though, and there are here.")

    rtree = rootTree(utree,outGroupTaxaL)
    
    Phylo.write(rtree, rtreeFN, 'newick')
