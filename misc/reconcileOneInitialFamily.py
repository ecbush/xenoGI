import sys,os
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import parameters,trees,Tree,genomes,families,Family,new_DTLOR_DP


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
    paramFN = sys.argv[2]
    ifamNum = int(sys.argv[3])
    D =  int(sys.argv[4])
    T =  int(sys.argv[5])
    L =  int(sys.argv[6])   
    O =  int(sys.argv[7])
    R =  int(sys.argv[8])

    # load trees
    speciesRtreeO = Tree.Rtree()
    speciesRtreeO.fromNewickFileLoadSpeciesTree(speciesTreeFN)

    # parameters file
    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)
    
    # ifam and support data
    genesO = genomes.genes(paramD['geneInfoFN'])
    initialFamiliesO = families.readFamilies(paramD['initFamilyFN'],speciesRtreeO,genesO,"initial")
    ifamO = initialFamiliesO.getFamily(ifamNum)
    tipMapD=families.getTipMapping(ifamO.geneTreeO,genesO)
    gtLocusMapD = getGtLocusMap(ifamO.geneTreeO,initialFamiliesO)

    # prepare input arguments, initFam number doesn't matter, so made it 1
    argT = (1,speciesRtreeO,ifamO.geneTreeO.unroot(),tipMapD,gtLocusMapD,D,T,L,O,R) 
    # reconcile
    initFamNum,geneRtreeO,graphD,minCost = families.reconcileOneUnRootedGeneTree(argT)

    # output
    
    print("Rooted gene tree:",geneRtreeO.toNewickStr())
    
    # put in initial family object so we can use its methods
    newIfamO = Family.initialFamily(None,None,geneTreeO=geneRtreeO,dtlorCost=minCost,dtlorGraphD=graphD)

    print()
    print("Num mprs",newIfamO.countMPRs())
    print()
    
    print("Printing all MPRs:")

    # We need to overwrite the DTLOR costs in paramD with the ones passed in by user
    paramD['duplicationCost'] = D
    paramD['transferCost'] = T
    paramD['lossCost'] = L
    paramD['originCost'] = O
    paramD['rearrangeCost'] = R
    for mprOrigFormatD,mprNodeFormatD in newIfamO.iterMprReconDFromGraph(speciesRtreeO.preorder(),paramD,False):
        
    
        print()
        print("Reconcilation")
        
        # put in origin family object so we can use its methods
        fam = Family.originFamily(None,None,geneTreeO=geneRtreeO,dtlorMprD=mprNodeFormatD)
        fam.printReconByGeneTree(genesO)
