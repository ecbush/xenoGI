import sys, os, statistics
sys.path.append(os.path.join(sys.path[0],'..'))
import parameters,genbank,blast,trees,genomes,Score,scores

# prints out average sizes of windows in bp, given the parameters
# we're using (which are in number of genes)

def calcWinSize(neighborT,geneNames,geneInfoD):
    '''Given a tuple of genes, find the size of the window they define in bp.'''
    # first one
    geneName = geneNames.numToName(neighborT[0])
    commonName,locusTag,descrip,chrom,start,end,strand = geneInfoD[geneName]
    minSt=int(start)
    maxEnd=int(end)

    # rest
    for geneNum in neighborT:
        geneName = geneNames.numToName(geneNum)
        commonName,locusTag,descrip,chrom,start,end,strand = geneInfoD[geneName]
        start=int(start)
        end=int(end)
        if start < minSt:
            minSt = start
        if end > maxEnd:
            maxEnd = end

    return maxEnd-minSt
        
        
def printWinSizeSummary(neighborTL):
    '''Given a list where index is genes and the values are neighbor genes, calculate the size of this window in bp for each gene. Return the mean and standard deviation.'''

    winL = []
    for neighborT in neighborTL:
        winL.append(calcWinSize(neighborT,geneNames,geneInfoD))

    median = statistics.median(winL)
    mean = statistics.mean(winL)
    stdev = statistics.stdev(winL)
                    
    print("  median",round(median))
    print("  mean",round(mean))
    print("  stdev",round(stdev))

## mods for core stuff (requires changing functions, so we move them here)

def createCoreSyntenyT(geneToAabrhT,geneOrderT,coreSynWsize):
    '''Create and return core synteny tuple. The index of this corresponds to gene number. The value at that index is a tuple of the aabrh numbers for the syntenic genes within coreSynWsize.'''

    coreSyntenyL =[None] * len(geneToAabrhT)
    geneAabrhOrderL = scores.createGeneAabrhOrderL(geneToAabrhT,geneOrderT)
    for element in geneAabrhOrderL:
        if element != None:
            # internal nodes are None, having no genes themselves
            for geneNumT,aabrhNumT in element:
                #print('h',len(geneNumT),len(aabrhNumT))
                for pos in range(len(geneNumT)):
                    geneNum=geneNumT[pos]
                    # changed to geneNumT here
                    coreSyntenyL[geneNum]=scores.getAabrhContext(pos,geneNumT,coreSynWsize)
    return tuple(coreSyntenyL)





if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    # here we can adjust if we want
    synWSize = paramD['synWSize']
    coreSynWsize = paramD['coreSynWsize']

    
    ## load data structures we'll use below
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])

    # an object for gene name conversions
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    geneOrderT=genomes.createGeneOrderTs(paramD['geneOrderFN'],geneNames,subtreeL,strainStr2NumD)

    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])

    aabrhL = scores.loadOrthos(paramD['aabrhFN'])
    
    #### synSc

    print("Stats for synSc windows.")
    neighborTL = scores.createNeighborL(geneNames,geneOrderT,synWSize)
    printWinSizeSummary(neighborTL)



    #### coreSynSc

    print("Stats for coreSynSc windows.")
    geneToAabrhT = scores.createGeneToAabrhT(aabrhL,geneNames)
    coreSyntenyT = createCoreSyntenyT(geneToAabrhT,geneOrderT,coreSynWsize)

    printWinSizeSummary(coreSyntenyT)
