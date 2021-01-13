import sys,os,glob
sys.path.insert(0,os.path.join(sys.path[0],'..'))
from xenoGI import xenoGI,trees,blast,scores,parameters,genomes
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
    scaffoldTree = trees.Rtree()
    scaffoldTree.fromNewickFileLoadSpeciesTree(paramD['scaffoldTreeFN'])
    outGroupL = [paramD['outGroup']]
    rootFocalClade = "".join([x for x in scaffoldTree.children(scaffoldTree.rootNode) if x != outGroupL[0]]) # get child that is not outgroup
    # because xenoGI.printAnalysisWrapper is also used by regular
    # xenoGI, we pass in tree file rather than tree here
    xenoGI.printAnalysisWrapper(paramD,paramD['scaffoldTreeFN'],rootFocalClade)
    xlMode.xlAnalysisSummary(paramD)

def interactiveAnalysisWrapper(paramD):
    """Enter interactive mode."""
    import code,sys
    
    ## funcs
    
    def geneToLocFam(geneNum):
        '''Given a gene number, return the locus family number that it was
mapped to when genes were mapped to the scaffold.'''
        lfNum = locFamNumAr[geneNum] 
        if lfNum == unmappedVal:
            return None
        else: return lfNum
        
    ## interactive main
    
    # read array
    locFamNumAr,unmappedVal = readAllGenesMapToScaffoldAr(paramD['allGenesMapToScaffoldFN'])
    
    # set up interactive console
    vars = globals()
    vars.update(locals())
    readline.set_completer(rlcompleter.Completer(vars).complete)
    readline.parse_and_bind("tab: complete")
    code.InteractiveConsole(vars).interact()
    
def debugWrapper(paramD):
    '''Take us into interactive mode for debugging.'''
    
    import code,sys
    from xenoGI import trees,genomes,families,islands,analysis
    
    # set up interactive console
    vars = globals()
    vars.update(locals())
    readline.set_completer(rlcompleter.Completer(vars).complete)
    readline.parse_and_bind("tab: complete")
    code.InteractiveConsole(vars).interact()

    
if __name__ == "__main__":

    ## check command line
    try:
        assert(len(sys.argv) == 3)
        paramFN=sys.argv[1]
        task = sys.argv[2]
        assert(task in ['parseGenbank', 'obtainCoreOrthoSets','makeSpeciesTree', 'makeScaffold', 'printAnalysisXL', 'runAll','interactiveAnalysis','debug'])
    
    except:
        print(
            """
   Exactly two arguments required.
      1. A path to a parameter file.
      2. The task to be run which must be one of: parseGenbank, obtainCoreOrthoSets, makeSpeciesTree, makeScaffold, printAnalysisXL, runAll or interactiveAnalysis.

   For example: 
      python3 runXlMode.py xlParams.py parseGenbank
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
        print("makeScaffold",file=sys.stderr)
        makeScaffoldWrapper(paramD)
        print("printAnalysisXL",file=sys.stderr)
        printAnalysisXLWrapper(paramD)

    #### interactiveAnalysis
    elif task == 'interactiveAnalysis':
        interactiveAnalysisWrapper(paramD)
        
    #### debug
    elif task == 'debug':
        debugWrapper(paramD)
