import sys,os
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,genbank,trees,genomes,scores,families,islands

## Does island making without family making (i.e. the second half of
## what xenoGI.py does).

if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    ## load data structures we'll use below
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])

    # an object for gene name conversions
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    ## read scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)

    ## load gene families
    familyL = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)
    
    ## group gene families into islands
    outputSummaryF = open(paramD['outputSummaryFN'],'w')
    islands.makeIslands(geneOrderT,geneNames,subtreeL,tree,paramD['proxThreshL'],familyL,paramD['numThreads'],strainStr2NumD,strainNum2StrD,paramD['rootFocalClade'],paramD['islandOutFN'],outputSummaryF)
    outputSummaryF.close()
