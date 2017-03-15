## This parameter file contains python expressions with info on files
## and parameters

#### Input files ####

# Tree file in newick format. This should have named internal
# nodes. It does not need to have branch lengths (if it has them, they
# will be ignored).
treeFN='testE.tre'

# unix style file path to genbank gbff files
genbankFilePath = 'ncbi/*.gbff'

#### Parse output files ####

# gene order
geneOrderFN = 'geneOrder.txt'

# list of protein names found to be redundant in genbank files (will not
# be used)
redundProtsFN = 'redundProts.txt'

# gene descriptions (for use in analysis)
geneInfoFN = 'geneInfo.txt'

# A file specifying the mapping between genbank file names and human
# readable names. These human readable names are then used to to refer
# to the various species in subsequent analysis. They should match
# what's in the tree. If we don't want to rename and will use NCBI
# file names in later analysis, set to None.
fileNameMapFN = 'ncbiHumanMap.txt'

# unix style file path to fasta files
fastaFilePath = 'fasta/*.fa'


#### Blast output ####

# blast command line (except for db,query and outfiles)
blastCLine = 'blastp -matrix BLOSUM62 -gapopen 11 -gapextend 1 -evalue 0.01 -seg yes -outfmt 6'

# unix style file path to blast output files
blastFilePath = 'blast/*.out'


#### Algorithm output files ####

# alignment based similarity scores file
rawScoresFN = 'rawScores.out'

# sets of all around best reciprocal hits
aabrhFN = 'aabrh.out'

# normalized scores file
normScoresFN = 'normScores.out'

# synteny based scores file
synScoresFN = 'synScores.out'

# family file
familyFN='fam.out'

# group file
groupOutFN = 'groups.out'

# file with summary info about family and group formation
outputSummaryFN = 'outputSummary.txt'

#### Algorithm parameters ####

# in parallel code, how many threads to use
numThreads = 50


# in calculation of normalized scores, we get set of all around best
# reciprocal hits. This is the evalue threshold used there.
evalueThresh = 0.001

# alignment parameters for making scores
# note, since parasail doesn't charge extend on the first base its
# like emboss, and not like ncbi blast. NCBI blast uses existance 11,
# extend 1 for blosum 62, thus we should use open 12 gapExtend 1 here.
gapOpen = 12
gapExtend = 1
matrix = 'parasail.blosum62'

# Synteny window size, that is the size of the neighborhood of each
# gene to consider when calculating synteny scores. (measured in
# number of genes)
synWSize = 30

# When calculating synteny score between two genes, the number of
# pairs of scores to take (and average) from the neighborhoods of
# those two genes
numSynToTake = 25

# Minimum normalized score for family formation. This should be used
# as an extreme lower bound, to eliminate those things that are so
# obiously dissimilar that they could not be homologous by descent
# from the node under consideration. In units of standard deviation,
# centered around 0.
minNormThresh = -7.0

# Minimum synteny score value we'll accept when putting a gene in a
# family. Also applies to seeds. Note synteny scores are based on
# normalized scores
minSynThresh = -2.0

# Synteny score threshold for using synteny to adjust a raw
# score. Setting this lower makes us use synteny more, and thus will
# tend to make us put more genes in families. This is a normScore type
# score
synAdjustThresh = 0

# We use syntenty scores to adjust similarity scores in family
# finding. This parameter specifies the amount we multiply a rawScore
# by during this adjustment.
synAdjustExtent = 1.05

# In deciding whether to merge two groups, we judge partly based on
# the proximity of their genes. The elements of this list are tuples
# of the form (proximity threshold, rscore level). For example (1,0)
# says we'll consider proximity to mean adjacency (proximity 1 means
# adjacent genes), and we'll join groups if the rscore values is 0 or
# above. The algorithm loops through the list, using the criteria in
# the first tuple first, then proceeding to the second if there is on
# and so on.
proxThreshL = [(1,0),(2,2)]


#### Visualization and analysis ####

# Creating gff files of groups for visualization in the IGB browser.
# We want to display different groups with different colors. We do
# this by giving different 'score' values to different groups. Scores
# can range from 1 to 1000, but we only give certain discrete scores.

scoreNodeMapD = {'i9':1, 'i8':10} # groups with these mrca values always get this score

# The score for rest is based on the list below. It was made with
# createGroupGffs.createPotentialScoresL(100,1001,200,50)
potentialScoresL=[100, 300, 500, 700, 900, 150, 350, 550, 750, 950, 200, 400, 600, 800, 1000, 250, 450, 650, 850]

#### Visualization and analysis output files ####

# unix style file path to gff output files
gffFilePath = 'gff/*-testE-group.gff'


