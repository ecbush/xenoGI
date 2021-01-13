## This parameter file contains python expressions with info on files
## and parameters.

## Note that many other parameters are specified in the parameters.py
## file in the repository. It possible to change the values for those
## other parameters simply by adding them (with a different value) to
## this file. Values specified here will override whatever is in
## parameters.py.

#########################
#### Basic xenoGI #######
#########################

#### Input Assembly files ####

# unix style file path to genbank gbff files
genbankFilePath = 'ncbi/*.gbff'

# A file specifying the mapping between genbank file names and human
# readable names. These human readable names are then used to to refer
# to the various species in subsequent analysis. They should match
# what's in the tree. If we don't want to rename and will use NCBI
# file names in later analysis, set to None.
fileNameMapFN = None

#### Blast ####

# absolute path to the directory containing the blastp and makeblastdb
# executables. On Windows you may need to put in a second slash as an
# escape, e.g. 'C:\\Users\\guest\\blast-2.7.1+\\bin'
blastExecutDirPath = '/usr/bin/'

#### Making species trees ####

# For the case if you don't have a species tree to begin with and
# intend to make one.

# full path to the ASTRAL executable. On Windows you may need to put
# in a second slash as an escape.
# e.g. 'C:\\Users\\guest\\Astral\\astral.5.6.3.jar'
astralPath = '/usr/local/Astral/astral.5.6.3.jar'

# Single outgroup species to be used in rooting species
# tree. Uncomment and enter value here if making species tree.

outGroup = 'GCF_000236925.1_ASM23692v1_genomic'

#### Making gene trees ####

# If we should use DNA based alignments to make tree, then this should
# be True, otherwise False
dnaBasedGeneTrees = True

# full paths to muscle, FastTree, java and ASTRAL. On Windows you may
# need to put in a second slash as an escape,
# e.g. 'C:\\Users\\guest\\muscle\\muscle.exe'
musclePath = '/usr/bin/muscle'
fastTreePath = '/usr/local/bin/FastTree'
javaPath = '/usr/bin/java'

#### Family formation ####
# DTLOR parameters should be integers. duplicationCost <= transferCost
duplicationCost = 2
transferCost = 6
lossCost = 1
originCost = 6
rearrangeCost = 7

# reconciliation with costs permissive to origin events
# reconcilePermissiveOriginGeneListPath specifies a file giving a list
# of xenoGI genes (one per line, string form) which we should use
# pemissive-origin reconciliation on. For each of these genes, we
# identify the initial family it belongs to, and then do
# reconciliation with permissive costs. By default, we set to None
# (and don't do this type of reconciliation). Users can override by
# uncommenting the line below, and creating the corresponding file
#reconcilePermissiveOriginGeneListPath = 'permissiveOriginGeneList.txt'


#### Parallelization ####

# in parallel code, how many threads to use
numProcesses = 50

#########################
#### xlMode specific ####
#########################

#### obtainCoreOrthoSets ####

# The number of randomly chosen genomes we use for finding the initial
# set of all around best reciprocal hit core genes
numGenomesInRandomSample = 15

# file to record all genomes in our random sample
genomesInRandomSampleFN = 'genomesInRandomSample.txt'

# sets of all around best reciprocal hits for initial random sample
randomSampleAabrhFN = 'randomSampleAabrh.out'

# name and location for fasta with sequence from all aabrh genes
randomSampleAabrhFastaFN = 'randomSampleAabrh.fa'

# set of core orthologs for all strains
allStrainCoreOrthosFN = 'allStrainCore.out'

#### makeSpeciesTree ####

# location for astral output
astralTreeFN = 'astral.tre'

# the file name for the final tree
speciesTreeFN = 'allStrains.tre'

#### Trim tree ####

# trimLeafNum specifies the number of strains in the scaffold tree
trimLeafNum = 15
scaffoldTreeFN = 'scaffold.tre'


#### Refine and map to scaffold ####

# User can specify what goes in the scaffold by pointing this to a
# file (one strain per line in file). This number of strains in this
# file should be less than trimLeafNum
userSpecifiedStrainsFileName = None

# database where we put the genes we chose to represent each family
scaffoldFamilyRepGenesFastaFN = 'scaffoldFamilyRepGenes.fa'
scaffoldFamilyRepGenesNumPerFamily = 4

# xlMapAlignCoverThresh is a threshold for the length of the blast alignment
# relative to query and subject length (ranges between 0 and 1)
xlMapAlignCoverThresh = 0.75

# location of binary file holding mapping of all genes onto scaffold families
allGenesMapToScaffoldFN = 'allGenesMapToScaffoldFam.bout'

# databases where we put all unmapped genes after the initial scaffold
numUnmappedGenesPerFasta = 3000
unMappedGenesFilePathStem = 'unMappedGenes'

# Number of strains to add to second scaffold
numStrainsToAddToScaffold = 5

#### Visualization and analysis output files ####

# xlMode specific analysis
xlAnalysisSummaryStem = 'xlAnalysisSummary'
