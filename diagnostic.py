# functions for checking our code, setting parameters etc.
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages

def scoreHists(scoreFN,outFN,numBins):
    '''Read through a scores file, and separate into all pairwise comparisons. Then plot hist of each.'''

    # currently, this seems to require a display for interactive
    # plots. would be nice to make it run without that...

    pairD = readScorePairs(scoreFN)

    pyplot.ioff() # turn off interactive mode
    with PdfPages(outFN) as pdf:
        for key in pairD:
            fig = pyplot.figure()
            pyplot.hist(pairD[key],bins=numBins)
            pyplot.title('-'.join(key))
            pdf.savefig()
            pyplot.close()

    #pyplot.show()
            

def readScorePairs(scoreFN):
    '''Read through a scores file, and separate into all pairwise
comparisons. Return as dict.'''
    
    pairD = {}

    f = open(scoreFN,'r')

    while True:
        s = f.readline()
        if s == '':
            break
        g1,g2,sc = s.rstrip().split('\t')
        sc = float(sc)

        sp1,restOfGene1 = g1.split('-')
        sp2,restOfGene2 = g2.split('-')

        key = tuple(sorted([sp1,sp2]))

        if key in pairD:
            pairD[key].append(sc)
        else:
            pairD[key] = [sc]
        
    f.close()
    return pairD
