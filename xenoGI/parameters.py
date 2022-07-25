# module for loading parameters.
import os

#### Base parameters

# These are parameters which users are not likely to need to change.
# We store them in a string of the same format as a parameters file.

baseParamStr = """

#### Genbank ####
# output from parsing genbank files

# gene order
geneOrderFN = 'geneOrder.txt'

# gene descriptions (for use in analysis)
geneInfoFN = 'geneInfo.txt'

# strain numbers and names
strainInfoFN = 'strainInfo.txt'

# unix style file path to fasta files
fastaFilePath = 'fasta/*.fa'

# listing of files with problems
problemGenbankFN = 'problemGenbankFiles.txt'


#### Blast ####

# unix style file path to blast output files
blastFilePath = 'blast/*.out'

# blast command line (except for value for evalue as well as db,query
# and outfiles)
blastCLine = 'blastp -matrix BLOSUM62 -gapopen 11 -gapextend 1 -seg yes -outfmt "6 qseqid sseqid evalue qlen qstart qend slen sstart send pident score" -evalue '

# Blast e-value threshold
evalueThresh = 1e-8

# alignCoverThresh is a threshold for the length of the blast alignment
# relative to query and subject length (ranges between 0 and 1)
alignCoverThresh = 0.65

# threshold for the percent identity in a hit
percIdentThresh = 0.35

# string to join the two strains compared in filename for blast
blastFileJoinStr = '_-VS-_'

#### Scores ####

# Note: for scores output files, if the extension we use here is
# .bout, the output will be saved in binary format. Otherwise it will
# be a less compact text based format
scoresFN = 'scores.bout'

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

# This helps determine the size of the arrays for our hash table of
# scores. It should be a value >= to 1. If one, it means the hash
# arrays have as many entries as there are score pairs.
hashArrayScaleFactor = 2

#### Making species trees ####

makeSpeciesTreeWorkingDir = 'makeSpeciesTreeWorkDir'
deleteSpeciesTreeWorkingDir = True # if True, we delete when done

# where to put gene trees for the aabrh hard core families
aabrhHardCoreGeneTreesFN = 'aabrhHardCoreGeneTrees.out'

# where to keep ASTRAL output. The final tree (which has been rooted,
# stripped of branch lengths and had internal nodes named) will be put
# in speciesTreeFN given above.
astralTreeFN = 'astralTree.tre'

#### Making gene trees ####

useGeneRaxToMakeSpeciesTrees = True # if False, use FastTree only

# directory to store trees for gene families
geneFamilyTreesDir = 'geneFamilyTrees'
deleteGeneFamilyTreesDir = True # if True, delete working dir when done

# file stems
aabrhHardCoreGeneTreeFileStem = 'aabrhHardCoreFam'
blastFamGeneTreeFileStem = 'blastFam'

## GeneRax

GeneRaxDNASubstModel = "GTR+G"
GeneRaxProtSubstModel = "LG+G"
GeneRaxReconcilationModel = "UndatedDTL"
GeneRaxSearchRadius = 5
GeneRaxMappingFileStem = "mapping_file_"
GeneRaxMappingFileExt = ".link"
GeneRaxInputListFN = "geneRaxInputList.txt"
GeneRaxOutputDirN = "geneRaxOutputDir"

#### Family formation ####

# Hard core all around best reciprocal hit gene families
aabrhFN = 'aabrhHardCore.out'

# Other families
blastFamilyFN = 'blastFam.out'
blastFamilyTreeFN = 'blastFamTrees.out'
initFamilyFN = 'initFam.out'
originFamilyFN = 'originFam.out'

# summary info about family formation
familyFormationSummaryFN = 'familyFormationSummary.out'

# Family formation involves several thresholds. We determine these by
# looking at histograms of scores.
scoreHistNumBins = 80

## Homology check
# parameters for finding the homologous (right) peak in raw score
# histograms. These put restrictions on the peaks we can
# find. widthRelHeight specifies where we measure the width of a peak,
# as a proportion of the way down from the top.
homologRightPeakLimit = 1.0
widthRelHeight = 0.9

# case 1 (normal case)
homologPeakWidthCase1 = 0.05
homologRequiredProminenceCase1 = 0.4
homologLeftPeakLimitCase1 = 0.65

# case 2 (extreme prominence. But allow to be very narrow)
homologPeakWidthCase2 = 0
homologRequiredProminenceCase2 = 6
homologLeftPeakLimitCase2 = 0.90

# case 3 (wide width with low prominence)
homologPeakWidthCase3 = 0.25
homologRequiredProminenceCase3 = 0.10
homologLeftPeakLimitCase3 = 0.65

# parameters for finding the non-homologous (left) peak in raw score
# histograms.
nonHomologPeakWidth = 0.15
nonHomologPeakProminence = 1
nonHomologLeftPeakLimit = 0
nonHomologRightPeakLimit = 0.6

## Splitting blast and initial families

# threshold of branch lengths for splitting unrooted gene trees from
# blast families. obtained by examining this distribution of maximum
# scores for aabrh hard core families.
quantileForObtainingSplitThresholds = .95
multiplierForObtainingSplitThresholds = 1.5

# The maximum size for blast and initial families is a multiple of the
# number of tips on the species tree. These parameters give the
# multiplier used. We want larger blast families than initial
# families, so the multiplier for blast should be larger.
maxBlastFamSizeMultiplier = 8
maxInitialFamSizeMultiplier = 4

# If some families still too large after splitting based on threshold,
# we force them to be split on a large internal branch. The choice of
# branch is based on branch length, and on achieving balance between
# the number of nodes on each side. This parameter affects how much
# weighting balance gets. Numbers > 1 mean we give more weight to
# balance.
forceSplitUtreeBalanceMultiplier = 10

## Synteny thresholds

# We determine synteny thresholds between pairs of strains based on
# the distribution of scores. quantileForObtainingSynThresholds gives
# the quantile we take from this dist. We then multiply this by
# multiplierForObtainingSynThresholds to get a threshold.
quantileForObtainingSynThresholds = 0.1
multiplierForObtainingSynThresholds = 0.5

#### LocusIsland formation ####

# island file
islandOutFN = 'islands.out'

# Summary info about island formation
islandFormationSummaryFN = 'islandFormationSummary.out'

# geneProximityRange tells us how far out we want to go from each gene
# in recording proximity for geneProximityD
geneProximityRange = 2

# In deciding whether to merge two islands, we judge partly based on
# the proximity of their genes. proximityThreshold defines what
# proximity is. A proximityThreshold of 1 means adjacent genes
proximityThresholdMerge = 1

# rscThreshold of 0 means that we merge islands if their rscore is 0
# or above
rscThresholdMerge = 0

# maxClusterSize is the maximum size that we make clusters of islands
# in the first step of the merging process
maxClusterSize = 50

#### Reconciliation ####

## reconciliation with costs permissive to origin events (for families
## that repeatedly insert in same syntenic location).

# File giving a list of xenoGI genes (one per line, string form) which
# we should use pemissive-origin reconciliation on. For each of these
# genes, we identify the initial family it belongs to, and then do
# reconciliation with permissive costs. By default, we set to None
# (and don't do this type of reconciliation). Users can override in
# the param.py file.
reconcilePermissiveOriginGeneListPath = None

# Costs to be used when we want to be permissive to origin events
DTLRcostPermissiveOrigin = 100
originCostPermissiveOrigin = 1

#### Family Refinement ####

# In family refinement, we consider alternate most parsimonious
# reconcilations. In some cases, the number of MPRs for a single
# family may be prohibitively large. upperNumMprThreshold specifies a
# limit beyond which we will randomly sample from the
# possibilities. This parameter also specifies the number of samples
# to be taken in that case.
upperNumMprThreshold = 20

# In order to test alternate MPRs, we must obtain nearby
# families. This parameter specifies how far to travel in each
# direction (in genes) in various genomes to collect those families.
geneProximityRangeRefineFamilies = 4

# In family refinement, we consider all origin families in islands
# this size or smaller
islandLenThresholdRefineFamilies = 2


#### Visualization and analysis ####

# analysis output
analysisDir = 'analysis'
islandsFNStem = 'islands'
genesFNstem = 'genes'

# output for browsers
bedFilePath = 'bed/*-island.bed' # unix style file path to bed output files

bedNumTries = 100 # number of random tries to find best coloring for islands

# The possible rgb values for bed files is based on the list below. 
potentialRgbL = ['245,130,48', '188,143,14','0,102,0','230,26,135','0,0,128', '145,30,180','0,255,255','128,0,0','0,255,0', '255,0,255','240,230,140','32,178,170','240,128,128','255,0,0','100,149,237','60,179,113','0,255,130','128,128,128','94,94,94','102,51,0']

"""


#### Functions
def createParametersD(baseParamStr,paramFN):
    '''Create and return a parameters dictionary. First parses the string
passed in as baseParamStr. This consists of parameters users are less
likely to modify. Then adds user specific parameters contained in the
file paramFN. Note that because the user parameters are put in the
parameters dictionary second, it is possible for a user to override
one of the base parameters simply by including that parameter in their
parameter file.

    '''
    paramD={}
    baseParamL = baseParamStr.split('\n')
    paramD = addParametersToD(baseParamL,paramD)

    f=open(paramFN,'r')
    userParamL = f.read().split('\n')
    paramD = addParametersToD(userParamL,paramD)
    f.close()
    
    return paramD
    
def addParametersToD(paramL,paramD):
    '''Given a list of lines (e.g. from a parameters file) add parameters
to paramD. Each line is a string. Some will be comments or blank
lines, and others will be python assignment statements. We use the
assignment statements to create entries in paramD.
    '''

    for s in paramL:
        if s == '' or s.lstrip() == '' or s.lstrip()[0]=='#':
            continue
        # it's not a blank line, a line with only whitespace or a
        # comment.  so it must be an assignment statement.
        key,value = s.rstrip().split('=')
        key = key.strip()
        value = value.strip()
        paramD[key] = eval(value)
    
    return paramD    
    
def loadFileNameMapD(fileNameMapFN,genbankFileList=None):
    '''Create a dictionary with mappings between genbank file names and
the human readable names we use in the tree. If fileNameMapFN contains
a string, we load that file and construct the mappings based on this.
Expects file with one species per line: genbank name + white space +
human name. If fileNameMapFN is None, we create the mappings between
the full file names (with path and extension) and the stem of the file
name, which in this case should correspond to what is in the input
tree.
    '''
    fileNameMapD = {}
    strainNamesL = []
    if fileNameMapFN == None:
        for fullPathFN in genbankFileList:
            fn = os.path.split(fullPathFN)[-1]
            stem = os.path.splitext(fn)[0]
            fileNameMapD[fn] = stem
            strainNamesL.append(stem)
    else:
        f = open(fileNameMapFN,'r')
        while True:
            s = f.readline()
            if s == '':
                break
            elif s[0].isspace():
                continue
            genbankStem,human = s.rstrip().split()
            fileNameMapD[genbankStem] = human
            strainNamesL.append(human)
    return fileNameMapD,tuple(strainNamesL)
