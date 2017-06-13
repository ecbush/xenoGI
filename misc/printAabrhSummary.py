import sys,os
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,trees,genomes,scores,analysis




if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    
    # get familyL etc.
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    # scores
    scoresO = scores.readScores(paramD['scoresFN'],geneNames)


    aabrhL = scores.loadOrthos(paramD['aabrhFN'])
    strainNamesL=sorted([strainNum2StrD[leaf] for leaf in trees.leafList(tree)])
    aabrhRawScoreSummmaryD=scores.getAabrhRawScoreSummmaryD(strainNamesL,aabrhL,scoresO,geneNames)

    print("Mean and standard deviation of raw scores between aabrh orthologs for pairs of species.")
    
    rowL=[]
    rowL.append(['Species 1','Species 2','Mean','Standard dev'])
    rowL.append(['---------','---------','----','------------'])
    for keyT,valT in aabrhRawScoreSummmaryD.items():
        row=[]
        row.extend(keyT)
        row.append(format(valT[0],'.3f'))
        row.append(format(valT[1],'.3f'))
        rowL.append(row)

    analysis.printTable(rowL)
