## This parameter file contains python expressions with info on files
## and parameters

#### Input Assembly files ####

# unix style file path to genbank gbff files
genbankFilePath = 'ncbi/*.gbff'

## Parse output files

# gene order
geneOrderFN = 'geneOrder.txt'

# list of protein names found to be redundant in genbank files (will not
# be used)
redundProtsFN = 'redundProts.txt'

# gene descriptions (for use in analysis)
geneInfoFN = 'geneInfo.txt'

# strain numbers and names
strainInfoFN = 'strainInfo.txt'

# A file specifying the mapping between genbank file names and human
# readable names. These human readable names are then used to to refer
# to the various species in subsequent analysis. They should match
# what's in the tree. If we don't want to rename and will use NCBI
# file names in later analysis, set to None.
fileNameMapFN = 'ncbiHumanMap.txt'

# unix style file path to fasta files
fastaFilePath = 'fasta/*.fa'

#### Trees ####

# Species tree file in newick format. This should have named internal
# nodes. It does not need to have branch lengths (if it has them, they
# will be ignored).
# If makeSpeciesTree will be run (the tree is to be calculated) then
# this is where the input tree will be put
treeFN='example.tre'

# The node and the branch leading to it define the focal clade where
# islands will be reconstructed. Everything outside of this will be
# treated as ougroups (e.g. won't merge islands there).
rootFocalClade = 'i2'

## Parameters for optionally making a species tree

# The parameters below are for the case if you don't have a species
# tree to begin with and intend to make one.

# If we should use DNA based alignments to make tree, then this should
# be True, otherwise False
dnaBasedSpeciesTree = False

# full paths to muscle, FastTree, java and ASTRAL. On Windows you may need to
# put in a second slash as an escape,
# e.g. 'C:\\Users\\guest\\muscle\\muscle.exe'
musclePath = '/usr/bin/muscle'
fastTreePath = '/usr/local/bin/FastTree'
javaPath = '/usr/bin/java'
astralPath = '/usr/local/Astral/astral.5.6.3.jar'

makeSpeciesTreeWorkingDir = 'makeSpeciesTreeWorkDir'
deleteSpeciesTreeWorkingDir = True # if True, we delete when done

# where to put gene trees for the aabrh hard core families
aabrhHardCoreGeneTreesFN = 'aabrhHardCoreGeneTrees.out'

# where to keep ASTRAL output. The final tree (which has been rooted,
# stripped of branch lengths and had internal nodes named) will be put
# in treeFN given above.
astralTreeFN = 'astralTree.tre'

# for rooting the tree
outGroupTaxaL = ['S_bongori']

#### Blast ####

# absolute path to the directory containing the blastp and makeblastdb
# executables. On Windows you may need to put in a second slash as an
# escape, e.g. 'C:\\Users\\guest\\blast-2.7.1+\\bin'
blastExecutDirPath = '/usr/bin/'

# Blast e-value threshold
evalueThresh = 1e-8

# unix style file path to blast output files
blastFilePath = 'blast/*.out'


#### Algorithm output files ####

# Note: for scores output files, if the extension we use here is
# .bout, the output will be saved in binary format. Otherwise it will
# be a less compact text based format
scoresFN = 'scores.bout'

# Where to put all around best reciprocal hit gene families (ie the
# "hard core").
aabrhFN = 'aabrhHardCore.out'

# family file
familyFN='fam.out'

# island file
islandOutFN = 'islands.out'

# file with summary info about family and island formation
familyFormationSummaryFN = 'familyFormationSummary.out'
islandFormationSummaryFN = 'islandFormationSummary.out'

#### Algorithm parameters ####

# in parallel code, how many threads to use
numThreads = 50

# alignment parameters for making scores
# note, since parasail doesn't charge extend on the first base its
# like emboss, and not like ncbi blast. NCBI blast uses existance 11,
# extend 1 for blosum 62, thus we should use open 12 gapExtend 1 here.
gapOpen = 12
gapExtend = 1
matrix = 'parasail.blosum62'

# Synteny window size, that is the size of the neighborhood of each
# gene to consider when calculating synteny scores. (measured in
# number of genes). We go half this distance in either direction.
synWSize = 8

# When calculating synteny score between two genes, the number of
# pairs of scores to take (and average) from the neighborhoods of
# those two genes
numSynToTake = 3

# Core synteny window size, that is the number of core genes from the
# neighborhood of each gene to consider when calculating core synteny
# scores. (measured in number of genes). We go half this distance in
# each direction.
coreSynWsize = 20

# Threshold for the core gene synteny score required for family
# formation. (These scores range from 0 to 1). Note that if this is
# set to higher than 0.5, then there will need to be at least one core
# gene on both sides. This would mean that things at the end of a
# contig won't get added to families because they lack a common core
# gene on that side.
minCoreSynThresh = 0.5

#### Visualization and analysis output files ####

# output for browsers
bedFilePath = 'bed/*-island.bed' # unix style file path to bed output files

# analysis output
analysisFilePath = 'analysis/*.out'
islandsSummaryStem = 'islandsSummary'
genesFNstem = 'genes'
