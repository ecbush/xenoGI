# Anne Clark, Alex Eng, Bo Lee, Eliot Bush
import sys, parseFuncs

dbFN=sys.argv[1]
matrixFN = sys.argv[2]
outputFN = sys.argv[3]

# get list of DB files
F = open(dbFN, 'r')
DB_FILES = [x.rstrip() for x in F.readlines()]
F.close()



def getPossibleGenes(rHits):
    '''return dict containing an entry for each gene in species 0 that
    has a best reciprocal hit in all other species. 
    format of an entry is gene:((sp#,gene),(sp#,gene),...)).
    species are in order in the tuple'''

    # this eliminates impossible orthologs becuase every gene that
    # has an ortholog in all species should appear in first species
    possibleGenes = {}
    firstRow = rHits[0]
    firstDict = firstRow[1]
    for gene in firstDict:
        isPossible = True
        hits = []
        for colInd in xrange(2, len(firstRow)):
            compareDict = firstRow[colInd]
            hit = compareDict.get(gene)
            if hit == None:
                isPossible = False
                break
            else:
                hits.append((colInd, hit))
        if isPossible:
            possibleGenes[gene] = [(1, firstDict[gene])] + hits
    return possibleGenes


def getOrthologs(possibleGenes, rHits):

    # final list of pairwise orthologous genes in format 
    # [[(sp#, gene), (sp#, gene),...],...]
    orth = []
    
    for gene in possibleGenes:
        hits = possibleGenes[gene]

        # check that each gene in the list of hits is best reciprocal 
        # hit of all following genes
        isOrth = True
        for hit1Index in xrange(len(hits)):
            hit1 = hits[hit1Index][1]

            # check that current gene in list is best reciprocal hit
            # of all following genes in list
            for hit2Index in xrange(hit1Index+1, len(hits)): # j always > i
                hit2 = hits[hit2Index][1]
                sp1sp2Dict = rHits[hit1Index+1][hit2Index+1]
                match = sp1sp2Dict.get(hit1)
                if match != hit2: # oh noes! not reciprocal!
                    isOrth = False
                    break

            # if any pair not reciprocal best hits, then we fail
            if not isOrth:
                break

        # if we get through all genes without finding that some aren't
        # reciprocal best hits, then we keep this list of genes
        if isOrth:
            orth.append([(0,gene)]+hits)
    return orth

def run():
    f = open(matrixFN, "r")
    # get matrix of pairwise reciprocal best hits
    rHits = parseFuncs.matrixFileReader(f)
    f.close()

    # get sets of genes in first species has a reciprocal best
    # hit in each other species
    possibleGenes = getPossibleGenes(rHits)

    # then, out of the possible genes, get those that are all
    # pairwise reciprocal bets hits
    orth = getOrthologs(possibleGenes, rHits)

    # create string of tab-delimited genes, with each set of orthologs
    # on a new line
    s = ""
    for geneL in orth:
        for spL in geneL:
            s += spL[1] + '\t'
        s += '\n'
    s=s[:-1] # remove last newline, so as not to have blank line at end
    return s


fout=open(outputFN,"w")

# print all of the species filenames
for i in xrange(len(DB_FILES)):
    print >>fout, i, DB_FILES[i]

# and then print all the lines of tab-delimited orthologous genes
print >>fout, run()

fout.close()
