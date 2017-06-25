
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

delProb = 0.8
minDel = 1
maxDel = 10

dupProb = 0.2
minDup = 1
maxDup = 10

invProb = 0.1
minInv = 100
maxInv = 500

hgtProb = 0.6
minHgt = 1
maxHgt = 10


## Sequence level evolution
# note the amount of evolution here really depends on the branch
# lengths in the tree, and doesn't depend on the step size (ie its
# different than all the genome scale events and small indels in that
# way).

model = "WAG" # used in pyvolve Markov model

smallIndelProb = 0 #0.003 # per gene per step
minIndel = 1 # min and max indel size in aa
maxIndel = 5

## Output files

logFile = 'logSim.txt'
genomeFilePrefix = 'sp'

# unix style file path to fasta files
fastaFilePath = 'fasta/*.fa'
