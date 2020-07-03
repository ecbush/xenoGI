import sys,os,copy
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import parameters,trees,genomes,families,Family


def getGtLocusMap(geneTree,initialFamiliesO):
    '''Get the mapping of genes to syntenic locations for this gene
tree. We rely on xenoGI having been run already to the family
formation state (to obtain an initialFamiliesO object).'''

    gtLocusMapD = {}
    for geneNum in trees.leafList(geneTree):
        gtLocusMapD[geneNum] = None

    # now fill with syntenic location numbers
    for lfO in initialFamiliesO.iterLocusFamilies():
        for geneNum in lfO.iterGenes():
            if geneNum in gtLocusMapD:
                gtLocusMapD[geneNum] = lfO.locusNum

    return gtLocusMapD

         
if __name__ == "__main__":
    
    speciesTreeFN = sys.argv[1]
    geneTreeFN = sys.argv[2]
    paramFN =  sys.argv[3]
    D =  float(sys.argv[4])
    T =  float(sys.argv[5])
    L =  float(sys.argv[6])   
    O =  float(sys.argv[7])
    R =  float(sys.argv[8])

    # load trees
    speciesTree = trees.readTree(speciesTreeFN)
    geneTree = trees.loadOneGeneTree(geneTreeFN) # expects xenoGI gene numbers on tips

    # parameters file
    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)
    
    # support data
    genesO = genomes.genes(paramD['geneInfoFN'])
    initialFamiliesO = families.readFamilies(paramD['initFamilyFN'],speciesTree,genesO)
    tipMapD=families.getTipMapping(geneTree,genesO)
    gtLocusMapD = getGtLocusMap(geneTree,initialFamiliesO)
    locusMapForRootingD = trees.createLocusMapForRootingD(geneTree,copy.deepcopy(gtLocusMapD))

    # prepare input arguments, initFam number doesn't matter, so made it 1
    argT = (1,speciesTree,geneTree,tipMapD,gtLocusMapD,locusMapForRootingD,D,T,L,O,R)
 
    # reconcile
    initFamNum,optRootedGeneTree,optMPR,minCost = families.reconcile(argT)

    # output
    reconD = families.convertReconBranchToNode(optMPR,optRootedGeneTree)
    
    print("Rooted gene tree:")
    newickOutTree = trees.tupleTree2NoBrLenNewick(optRootedGeneTree)
    print(newickOutTree)
    print()

    # put in family object so we can use its methods
    fam = Family.Family(None,None,optRootedGeneTree,reconD)
    print("Reconciliation:")
    fam.printReconByGeneTree()
