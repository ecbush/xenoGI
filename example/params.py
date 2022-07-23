## This parameter file contains python expressions with info on files
## and parameters.

## Note that many other parameters are specified in the parameters.py
## file in the repository. It possible to change the values for those
## other parameters simply by adding them (with a different value) to
## this file. Values specified here will override whatever is in
## parameters.py.

#### Input Assembly files ####

# unix style file path to genbank gbff files
genbankFilePath = 'ncbi/*.gbff'

# A file specifying the mapping between genbank file names and human
# readable names. These human readable names are then used to to refer
# to the various species in subsequent analysis. They should match
# what's in the tree. If we don't want to rename and will use NCBI
# file names in later analysis, set to None.
fileNameMapFN = 'ncbiHumanMap.txt'

#### Trees ####

# Species tree file in newick format. This should have named internal
# nodes. It does not need to have branch lengths (if it has them, they
# will be ignored).
# If makeSpeciesTree will be run (the tree is to be calculated) then
# this is where the input tree will be put
speciesTreeFN='example.tre'

# The node and the branch leading to it define the focal clade where
# islands will be reconstructed. Everything outside of this will be
# treated as ougroups (e.g. we won't merge islands there).
rootFocalClade = 's2'

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
astralPath = '/usr/local/bin/astral'

# Single outgroup species to be used in rooting species
# tree. Uncomment and enter value here if making species tree.

#outGroup = ''

#### Making gene trees ####

# If we should use DNA based alignments to make tree, then this should
# be True, otherwise False
dnaBasedGeneTrees = True

# full paths to muscle, FastTree, java and ASTRAL. On Windows you may
# need to put in a second slash as an escape,
# e.g. 'C:\\Users\\guest\\muscle\\muscle.exe'
musclePath = '/usr/local/bin/muscle'
fastTreePath = '/usr/local/bin/fasttree'
geneRaxPath = '/usr/local/bin/generax'
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
