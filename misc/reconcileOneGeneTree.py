import sys,os,copy
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import parameters,trees,Tree,genomes,families,Family


def getGtLocusMap(geneUtreeO,initialFamiliesO):
    '''Get the mapping of genes to syntenic locations for this gene
tree. We rely on xenoGI having been run already to the family
formation state (to obtain an initialFamiliesO object).'''

    gtLocusMapD = {}
    for geneNumStr in geneUtreeO.leaves():
        gtLocusMapD[geneNumStr] = None

    # now fill with syntenic location numbers
    for lfO in initialFamiliesO.iterLocusFamilies():
        for geneNum in lfO.iterGenes():
            if str(geneNum) in gtLocusMapD:
                gtLocusMapD[str(geneNum)] = lfO.locusNum

    return gtLocusMapD

         
if __name__ == "__main__":
    
    speciesTreeFN = sys.argv[1]
    geneTreeFN = sys.argv[2]
    paramFN =  sys.argv[3]
    D =  int(sys.argv[4])
    T =  int(sys.argv[5])
    L =  int(sys.argv[6])   
    O =  int(sys.argv[7])
    R =  int(sys.argv[8])

    # load trees
    speciesRtreeO = Tree.Rtree()
    speciesRtreeO.fromNewickFileLoadSpeciesTree(speciesTreeFN)

    geneUtreeO = Tree.Utree()
    geneUtreeO.fromNewickFile(geneTreeFN) # expects xenoGI gene numbers on tips

    # parameters file
    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)
    
    # support data
    genesO = genomes.genes(paramD['geneInfoFN'])
    initialFamiliesO = families.readFamilies(paramD['initFamilyFN'],speciesRtreeO,genesO)
    tipMapD=families.getTipMapping(geneUtreeO,genesO)
    gtLocusMapD = getGtLocusMap(geneUtreeO,initialFamiliesO)

    # prepare input arguments, initFam number doesn't matter, so made it 1
    argT = (1,speciesRtreeO,geneUtreeO,tipMapD,gtLocusMapD,D,T,L,O,R) 
    # reconcile
    initFamNum,optGeneRtreeO,optMPR,minCost = families.reconcile(argT)
    
    # output
    reconD = families.convertReconBranchToNode(optMPR)
    
    print("Rooted gene tree:")
    newickOutTree = optGeneRtreeO.toNewickStr()
    print(newickOutTree)
    print()

    # put in family object so we can use its methods
    fam = Family.Family(None,None,optGeneRtreeO,reconD)
    print("Reconciliation:")
    fam.printReconByGeneTree()
