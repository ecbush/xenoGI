
## tree

# based on testB.tre
tree = (14, (12, (10, (6, (4, (2, (0, (), (), 0.02), (1, (), (), 0.02), 0.04), (3, (), (), 0.06), 0.04), (5, (), (), 0.1), 0.1), (9, (7, (), (), 0.04), (8, (), (), 0.04), 0.16), 0.1), (11, (), (), 0.3), 0.1), (13, (), (), 0.4), 0)

## Genome basics

geneNumber = 200 
minSeqLen = 200
maxSeqLen = 300 

step = .0025 # the amount of branch length in a single step. A step is
             # our discrete unit for various genome scale events

## Probs of verious genome scale events occurring in one step. Note
## the probability of an event does not depend on the number of genes.
## Note that we want del to cancel dup + hgt

delProb = 0.4
minDel = 1
maxDel = 10

dupProb = 0.2
minDup = 1
maxDup = 10

invProb = 0.2
minInv = 5
maxInv = 10

hgtProb = 0.2
minHgt = 1
maxHgt = 10


## Sequence level evolution
# note the amount of evolution here really depends on the branch
# lengths in the tree, and doesn't depend on the step size (ie its
# different than all the genome scale events and small indels in that
# way).

model = "WAG" # used in pyvolve Markov model

smallIndelProb = 0.005 # per gene per step
minIndel = 1 # min and max indel size in aa
maxIndel = 5

## Output files

logFile = 'logSim.txt'
genomeFilePrefix = 'species'
