## This parameter file contains python expressions with info on files
## and parameters

#### Input files ####

# tree file
treeFN='example.tre'

# unix style file path to genbank gbff files
genbankFilePath = 'genbank/*.gbff'

#### Parse output files ####

# gene order
geneOrderFN = 'geneOrder.txt'

# list of protein names found to be redundant in genbank files (will not
# be used)
redundProtsFN = 'redundProts.txt'

# gene descriptions (for use in analysis)
geneDescriptionsFN = 'geneDescriptions.txt'

# unix style file path to fasta files
fastaFilePath = 'fasta/*.fa'


#### Blast output ####

#blastp -matrix BLOSUM62 -gapopen 11 -gapextend 1 -evalue 0.01 -seg yes -outfmt 6 -db fasta/Citrobacter.fa -query fasta/Klebsiella.fa -out aa

# blast command line (except for db,query and outfiles)
blastCLine = 'blastp -matrix BLOSUM62 -gapopen 11 -gapextend 1 -evalue 0.01 -seg yes -outfmt 6'

# unix style file path to blast output files
blastFilePath = 'blast/*.out'



#### Algorithm output files ####

# alignment based similarity scores file
scoresFN = 'scores.out'

# synteny based scores file
synScoresFN = 'synScores.out'

# family file
familyFN='fam.out'

# group file
groupOutFN = 'groups.out'




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

# We use syntenty scores to adjust similarity scores in family
# finding. This parameter specifies the maximum amount we can
# add to sim scores
synAdjustMaxExtent = 0.1

# Threshold for group scores. Below this we do not merge groups. These
# scores theoretically range from -2 to 2 inclusive.
groupScoreThreshold = 0

# alignment parameters for making scores
# note, since parasail doesn't charge extend on the first base its
# like emboss, and not like ncbi blast. NCBI blast uses existance 11,
# extend 1 for blosum 62, thus we should use open 12 gapExtend 1 here.
gapOpen = 12
gapExtend = 1
matrix = 'parasail.blosum62'
