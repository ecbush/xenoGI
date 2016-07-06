# Anne Clark, Alex Eng, Bo Lee (??), Eliot Bush
import sys

NO_THRESHOLD = False

dbFN=sys.argv[1]
outputFN = sys.argv[2]
blastDir = sys.argv[3] # path to blast files
IDENTITY_THRESHOLD = int(sys.argv[4])
LENGTH_THRESHOLD = float(sys.argv[5]) # lengths should be within this
                                      # percent of each other

# get list of strain names
F = open(dbFN, 'r')
STRAIN_NAMES = [x.rstrip() for x in F.readlines()]
F.close()

# Usage:
# python ebcode/reciprocalMatrix.py full/dbList.txt full/recipMatrix.out blast/out/ 85 .2


def getReciprocalHits(strainName1, strainName2):
    """a file contains the results of blasting one species
    against another; returns a two-element list that contains
    dictionaries for the first and second species; the key in each
    dictionary is a gene, and the value is that gene's best
    reciprocal hit from the other species; all keys in a dictionary
    are from the same species"""

    # get best hits of 1 blasted against 2...
    hits1 = getHits(blastDir+strainName1+'-'+strainName2+'.out')
    # ...and 2 blasted against 1
    hits2 = getHits(blastDir+strainName2+'-'+strainName1+'.out')

    # then store only the reciprocal best hits
    recHits1 = {}
    for query in hits1:
        hit = hits1[query]
        nextHit = hits2.get(hit)
        if nextHit == query:
            recHits1[query] = hit

    return recHits1

def getHits(fileName):
    """Given a BLAST output file, returns a dictionary keyed by the genes
    in the query species, with the values being the top hit (if any)
    for those genes. Assumes the blast hits for each query are given
    from most to least significant. which appears to be the
    case. Thresholds for minimum similarity and maximum length
    difference are globally defined.

    """

    f = open(fileName, 'r')
    hits = {}

    queryGene = ""

    # read through entire file
    line = f.readline()

    while line != '':
        
        # break the line up into its components
        L = line.split('\t')

        if len(L)<12:
            # its not an info line (likely header)
            line = f.readline()
            continue

        # gene names come first
        queryGene = L[0]
        hit = L[1]

        # percent id is the third element
        identity = int(float(L[2]))

        # lengths calculated from their start and end indices
        queryLength = int(L[7]) - int(L[6])
        hitLength = int(L[9]) - int(L[8])

        # store the query and hit if they meet thresholds
        if NO_THRESHOLD or (identity >= IDENTITY_THRESHOLD and \
                                similarLengths(queryLength, hitLength)):
            hits[queryGene] = hit

        # we only want the first hit for any query gene
        while line.split('\t')[0] == queryGene:
            line = f.readline()

    f.close()
    return hits

def similarLengths(length1, length2):
    """Are the lengths less than LENGTH_THRESHOLD% different?"""

    return (1 - (float(min(length1, length2)) / max(length1, length2))) < LENGTH_THRESHOLD

def getAllReciprocalHits():
    '''return an upper-diagonal (N-1)xN matrix where each entry
    [i][j] (j > i) contains a dictionary of best reciprocal hits
    between species i and species j, as indexed in list STRAIN_NAMES;
    key is always gene in species i, and value in species j'''

    rHits = []
    
    # create the matrix but skip last row because it will just have
    # Sp_n vs Sp_n
    for i in xrange(len(STRAIN_NAMES)): 

        rHits.append([])

        # fill in the rest of the row with the dictionary of best
        # reciprocal hits between species i and j (keyed by species i)
        for j in xrange(len(STRAIN_NAMES)):
            if j != i:
                rHits[i].append(getReciprocalHits(STRAIN_NAMES[i], STRAIN_NAMES[j]))
            else:
                rHits[i].append(None)
    return rHits

def run():
    f=open(outputFN,"w")
    rHits = getAllReciprocalHits()
    for x in range(len(rHits)):
        print >>f, '$$$'
        for y in range(len(rHits[x])):
            if rHits[x][y] != None:
                print >>f, ">" + STRAIN_NAMES[x] + "#%#"  + STRAIN_NAMES[y]
                for z in rHits[x][y].keys():
                    print >>f, z + "#%#" + rHits[x][y][z]
            else:
                print >>f, '!!!'
    f.close()

run()
