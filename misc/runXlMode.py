import sys, os, glob
sys.path.insert(0,os.path.join(sys.path[0],'../xenoGI'))
from xenoGI import xenoGI,trees, blast, scores, parameters, fasta, genomes
import xlMode

## Wrappers

def makeSpeciesTreeWrapper(paramD):
    '''Starting with a set of core orthologs for all strains, create a species tree.'''
    genesO = genomes.genes(paramD['geneInfoFN'])
    allStrainCoreL = scores.loadOrthos(paramD['allStrainCoreOrthosFN'])
    trees.makeSpeciesTree(paramD,allStrainCoreL,genesO)
    
def makeScaffoldWrapper(paramD):
    '''Starting with a species tree for all strains, obtain a
representative scaffold tree.
    '''
    xlMode.trimTree(paramD)
    xlMode.scaffoldWrapper(paramD)

def printAnalysisXLWrapper(paramD):
    '''Wrapper to run analysis.
    '''
    scaffoldTree = trees.readTree(paramD['scaffoldTreeFN'])
    rootFocalClade = trees.getRootFocalCladeFromOutgroup(scaffoldTree,paramD['outGroup'])
    # because xenoGI.printAnalysisWrapper is also used by regular
    # xenoGI, we pass in tree file rather than tree here
    xenoGI.printAnalysisWrapper(paramD,paramD['scaffoldTreeFN'],rootFocalClade)
    xlMode.xlAnalysisSummary(paramD)

def debugWrapper(paramD):
    '''Take us into interactive mode for debugging.'''
    
    import code,sys
    from xenoGI import trees,genomes,families,islands,analysis

    
    code.interact(local=locals())

    
if __name__ == "__main__":

    ## check command line
    try:
        assert(len(sys.argv) == 3)
        paramFN=sys.argv[1]
        task = sys.argv[2]
        assert(task in ['parseGenbank', 'obtainCoreOrthoSets','makeSpeciesTree', 'makeScaffold', 'printAnalysisXL', 'runAll', 'debug'])
    
    except:
        print(
            """
   Exactly two arguments required.
      1. A path to a parameter file.
      2. The task to be run which must be one of: parseGenbank, obtainCoreOrthoSets, makeSpeciesTree, makeScaffold, printAnalysisXL or runAll.

   For example: 
      python3.7 fastXenoGI.py params.py parseGenbank
"""
            ,file=sys.stderr)
        sys.exit(1)
        
    #### load parameters and some other data we'll use below
    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN) # read in parameters

    #### parseGenbank
    if task == 'parseGenbank':
        xenoGI.parseGenbankWrapper(paramD)

    #### obtainCoreOrthoSets
    elif task == 'obtainCoreOrthoSets':
        xlMode.obtainCoreOrthoSetsWrapper(paramD)

    #### makeSpeciesTree
    elif task == 'makeSpeciesTree':
        makeSpeciesTreeWrapper(paramD)
        
    #### makeScaffold
    elif task == 'makeScaffold':
        makeScaffoldWrapper(paramD)
        
    #### printAnalysisXL
    elif task == 'printAnalysisXL':
        printAnalysisXLWrapper(paramD)
        
    #### runAll
    elif task == 'runAll':
        xenoGI.parseGenbankWrapper(paramD)
        print("obtainCoreOrthoSets",file=sys.stderr)
        xlMode.obtainCoreOrthoSetsWrapper(paramD)
        print("makeSpeciesTree",file=sys.stderr)
        makeSpeciesTreeWrapper(paramD)
        print("makeRepTree",file=sys.stderr)
        makeScaffoldWrapper(paramD)
        print("printAnalysisXL",file=sys.stderr)
        printAnalysisXLWrapper(paramD)

    #### debug
    elif task == 'debug':
        debugWrapper(paramD)
