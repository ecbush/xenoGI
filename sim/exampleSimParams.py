
## tree

treeFN = 'testEBranch.tre'

## Genome basics

initialGeneNumber = 4000
minSeqLen = 157 # based on lengths in testE. 25 percentile
maxSeqLen = 402 # 75percentile

step = .0006 # the amount of branch length in a single step. A step is our discrete unit for various genome scale events

## Probs of verious genome scale events occurring in one step. Note
## the probability of an event does not depend on the number of genes.
## Note that we want del to cancel dup + hgt

delProb = 0.7
minDel = 1
maxDel = 50

dupProb = 0.35
minDup = 1
maxDup = 48

invProb = 0.1
minInv = 5
maxInv = 150

hgtProb = 0.35
minHgt = 2
maxHgt = 51


## Sequence level evolution
# note the amount of evolution here really depends on the branch
# lengths in the tree, and doesn't depend on the step size (ie its
# different than all the genome scale events and small indels in that
# way).

model = "WAG" # used in pyvolve Markov model

smallIndelProb = 0.00000015 # per gene per step
minIndel = 1 # min and max indel size in aa
maxIndel = 4

## Output files

logFile = 'logSim.txt'
genomeFilePrefix = 'sp'

# unix style file path to fasta files
fastaFilePath = 'fasta/*.fa'
