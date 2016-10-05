## This parameter file contains python expressions with info on files and parameters

#### Files ####

# tree file
treeFN='testA.tre'

# gene order
geneOrderFN = 'geneOrder.txt'

# alignment based similarity scores file
scoresFN = 'scores.out'

# synteny based scores file
synScoresFN = 'synScores.out'

# family file
familyFN='fam.out'

# group file
groupOutFN = 'groups.out'

# gene descriptions (for use in analysis)
geneDescriptionsFN = 'geneDescriptions.csv'


#### Parameters ####

# in parallel code, how many threads to use
numThreads = 50

# Synteny window size, that is the size of the neighborhood of each
# gene to consider when calculating synteny scores. (measured in
# number of genes)
synWSize = 20

# when calculating synteny score between two genes, the number of
# pairs of scores to take (and average) from the neighborhoods of
# those two genes
numSynToTake = 10

# minimum synteny score value we'll accept when putting a gene in a
# family. Also applies to seeds.
minSynThresh = 0.5

# Synteny score threshold for using synteny to adjust a similarity
# score. Setting this lower makes us use synteny more, and thus will
# tend to make us put more genes in families.
synAdjustThresh = 0.8

# We use a syntenty scores to adjust similarity scores in family
# finding. This parameter specifies the maximum amount we can
# add to sim scores
synAdjustMaxExtent = 0.1


# Threshold for group scores. Below this we do not merge groups. These
# scores theoretically range from -2 to 2 inclusive.
groupScoreThreshold = 0
