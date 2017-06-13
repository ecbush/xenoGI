import sys,os
sys.path.append(os.path.join(sys.path[0],'../'))
import parameters,genomes,trees,scores
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages



# Plots of scores

def scoreHists(scoresFN,outFN,numBins,geneNames,scoreType):
    '''Read through a scores file, and separate into all pairwise comparisons. Then plot hist of each.'''

    # currently, this seems to require a display for interactive
    # plots. would be nice to make it run without that...

    pairD = readScorePairs(scoresFN,geneNames,scoreType)

    pyplot.ioff() # turn off interactive mode
    with PdfPages(outFN) as pdf:
        for key in pairD:
            fig = pyplot.figure()
            pyplot.hist(pairD[key],bins=numBins)
            pyplot.title('-'.join(key))
            pdf.savefig()
            pyplot.close()


def readScorePairs(scoresFN,geneNames,scoreType):
    '''Read through a scores file, and separate into all pairwise
comparisons. Return as dict.'''
    
    pairD = {}

    scoresO = scores.readScores(scoresFN,geneNames=None)
    
    for gn1,gn2 in scoresO.iterateEdgesByEndNodes():
        sc = scoresO.getScoreByEndNodes(gn1,gn2,scoreType)
        sp1 = geneNames.numToStrainName(gn1)
        sp2 = geneNames.numToStrainName(gn2)
        key = tuple(sorted([sp1,sp2]))
        
        if key in pairD:
            pairD[key].append(sc)
        else:
            pairD[key] = [sc]
        
    return pairD


if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)
    scoreType = sys.argv[2]
    outFN = sys.argv[3]
    
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    scoresO = scores.readScores(paramD['scoresFN'],geneNames)

    scoreHists(paramD['scoresFN'],outFN,80,geneNames,scoreType)
