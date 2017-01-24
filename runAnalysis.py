import sys
import analysis,genomes,trees,families,scores,groups

if __name__ == "__main__":

    paramFN=sys.argv[1]

    params = __import__(paramFN.replace('.py', ''))

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)
    
    # load groups
    groupByNodeL=groups.readGroups(params.groupOutFN,tree,strainStr2NumD)
    
    # get familyT etc.
    geneNames = genomes.geneNames(params.geneOrderFN,strainStr2NumD,strainNum2StrD)

    
    geneInfoD = genomes.readGeneInfoD(params.geneInfoFN)

    familyT = families.readFamilies(params.familyFN,tree,geneNames,strainStr2NumD)

    
    gene2FamD=analysis.createGene2FamD(familyT)
    fam2GroupD=analysis.createFam2GroupD(groupByNodeL)

    # subtree list
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()


    geneOrderT=genomes.createGeneOrderTs(params.geneOrderFN,geneNames,subtreeL,strainStr2NumD)

    # scores
    rawScoresG = scores.readGraph(params.rawScoresFN,geneNames)
    normScoresG = scores.readGraph(params.normScoresFN,geneNames)
    synScoresG = scores.readGraph(params.synScoresFN,geneNames)
    
