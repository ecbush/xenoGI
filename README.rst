======
xenoGI
======

Code for reconstructing genome evolution in clades of microbes.

Requirements
------------

* NCBI blast+

  We need blastp and makeblastdb executables (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

* MUSCLE V5 (https://www.drive5.com/muscle/). For creating protein or DNA alignments.

* FastTree (http://www.microbesonline.org/fasttree/). For making gene trees.

* GeneRax (https://github.com/BenoitMorel/GeneRax). For making (species tree aware) gene trees. This is optional but recommended.
  
* Python 3

* Python package dependencies

  - Biopython (http://biopython.org/). This is for parsing genbank files and can be installed using pip:
      ``pip3 install biopython``

  - Parasail (https://github.com/jeffdaily/parasail). This is an optimized alignment library, used in calculating scores between proteins. It can also be installed using pip:
      ``pip3 install parasail``

  - Numpy (http://www.numpy.org/).
      ``pip3 install numpy``
    
  - Scipy (https://www.scipy.org/).
    ``pip3 install scipy``

(The pip you use needs to correspond to a version of Python 3. In some cases it may just be called pip instead of pip3).

* Additional dependencies

  If you make use of the ``makeSpeciesTree`` flag or ``xlMode.py``, you will also need the following

  - ASTRAL (https://github.com/smirarab/ASTRAL/).

* Comments on platforms.

  xenoGI is developed on Linux. The docker image (linked below) is the easiest way to run on Mac and Windows.

Installation
------------

Via pip::

  pip3 install xenoGI

(You will separately need to install blast+, MUSCLE, FastTree, and optionally GeneRax and ASTRAL.)

Via docker. For some instructions on using docker, go here:

  https://hub.docker.com/r/ecbush/xenogi

Using docker, xenoGI get's run within a virtual machine. This is nice because you don't have to worry about all the dependencies above (they're provided in our image). This does come at some cost in terms of performance.
  
Citation
--------

If you use xenoGI in a publication, please cite the following:

Bush EC, Clark AE, DeRanek CA, Eng A, Forman J, Heath K, Lee AB, Stoebel DM, Wang Z, Wilber M, Wu H. xenoGI: reconstructing the history of genomic island insertions in clades of closely related bacteria. BMC Bioinformatics. 19(32). 2018.

Liu J, Mawhorter R, Liu I, Santichaivekin S, Bush E, Libeskind-Hadas R. Maximum Parsimony Reconciliation in the DTLOR Model. BMC Bioinformatics. 22(394). 2021.

How to use
----------

An ``example/`` directory is included in this repository.

The sections below give some instructions about how to run xenoGI on this example. You can use this to make sure you've installed it properly and so forth. The github repository also contains a TUTORIAL which you can run through after completing the README.

The basic method works on a set of species with known phylogenetic relationships. In the example, these species are: E. coli K12, E. coli ATCC 11775, E. fergusonii and S. bongori. In cases where you don't know the species tree, xenoGI has methods to help you reconstruct it.

Required files
~~~~~~~~~~~~~~

The working directory must contain:

* A parameter file. In the provided ``example/`` directory this is called ``params.py``.

* A newick format tree representing the relationships of the strains. In the example this is called ``example.tre``. Note that branch lengths are not used in xenoGI, and ``example.tre`` does not contain branch lengths. Also note that internal nodes should be given names in this tree. In the example.tre we label them s0, s1 etc. The parameter ``speciesTreeFN`` in ``params.py`` has the path to this tree file. If a strain tree is not available, xenoGI has some accessory methods, described below, to help obtain one.

* A subdirectory of sequence files. In the example, this is called ``ncbi/``. Contained in this subdirectory will be genbank (gbff) files for the species. The parameter ``genbankFilePath`` in ``params.py`` has the path to these files.

Naming of strains
~~~~~~~~~~~~~~~~~

The system needs a way to connect the sequence files to the names used in the tree.

In the example, the sequence files have names corresponding to their assembly accession number from ncbi. We connect these to the human readable names in example.tre using a mapping given in the file ``ncbiHumanMap.txt``. This file has two columns, the first giving the name of the genbank file, and the second giving the name for the strain used in the tree file. In ``params.py`` the parameter ``fileNameMapFN`` is set to point to this file.

Note that the strain names should not contain any dashes, spaces, commas or special characters.

Another approach is to change the names of the sequence files to match what's in the tree. If you do this, then you should set ``fileNameMapFN = None`` in ``params.py``. (This is not necessary in the example, which is already set to run the other way).

Pointing xenoGI to various executables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before running xenoGI you'll have to ensure that it knows where various executables are. Edit ``params.py`` using a text editor such as emacs, vim, nano, Visual Studio Code etc. You should edit the following to give the absolute (full) path to the directory where the ``blastp`` and ``makeblastdb`` executables reside::

  blastExecutDirPath = '/usr/bin/'

(Change '/usr/bin/' to correspond to the right location on your system).

Also make sure that the absolute paths to MUSCLE and FastTree are correct in ``params.py`` (the parameters ``musclePath`` and ``fastTreePath``). If you intend to use generax to make species tree aware gene trees, then you also need to set ``geneRaxPath``. (The default parameter file is set to use generax, so unless you change the ``useGeneRaxToMakeSpeciesTrees`` parameter, described below, you'll need to supply a ``geneRaxPath``).

If you will be using the makeSpeciesTree functionality, then you will also need to specify ``astralPath`` and ``javaPath``.

Running the code
~~~~~~~~~~~~~~~~

If you install via pip, then you should have an executable script in your path called xenoGI.

You run the code from within the working directory. To run the example, you would cd into the ``example/`` directory. You will need to ensure that the ``params.py`` parameters file contains the correct path to the directory with the blastp and makeblastdb executables in it, as well as the MUSCLE and FastTree executables. Then, the various steps of xenoGI can be run all at once like this::

  xenoGI params.py runAll

They can also be run individually::

  xenoGI params.py parseGenbank
  xenoGI params.py runBlast
  xenoGI params.py calcScores
  xenoGI params.py makeFamilies
  xenoGI params.py makeIslands
  xenoGI params.py refine
  xenoGI params.py printAnalysis
  xenoGI params.py createIslandBed

If for some reason you don't want to install via pip, then you can download the repository and run the code like this::

  python3 path-to-xenoGI-github-repository/xenoGI-runner.py params.py runAll

(In this case you will have to make sure all the python package dependencies are satisfied.)

What the steps do
~~~~~~~~~~~~~~~~~

* ``parseGenbank`` runs through the genbank files and produces input files that are used by subsequent code. This step pulls out every CDS feature that has a ``/translation`` tag. The fields that are recorded (if present) are locus_tag, protein_id, product (that is gene description), and chromosomal coordinates as well as the protein sequence. If the parameter ``dnaBasedGeneTrees`` is True, the DNA sequence for each gene is kept as well.
  
* ``runBlast`` does an all vs. all protein blast of the genes in these strains. The number of processes it will run in parallel is specified by the ``numProcesses`` parameter in the parameter file. Before running a particular comparison, runBlast checks to see if the output file for that comparison already exists (e.g. from a previous run). If so it skips the comparison.
  
* ``calcScores`` calculates similarity and synteny scores between genes in the strains. It is also (mostly) parallelized.
  
* ``makeFamilies`` calculates gene families using blast, FastTree, GeneRax (optionally), and a customized variant of the DTL reconciliation algorithm called DTLOR. This approach considers synteny in the family formation process.

* ``makeIslands`` groups families according to their origin, putting families with a common origin together as islands. It is partly parallelized.

* ``refine`` reconsiders certain families in light of the output of makeIslands. In particular, this step looks at cases where there are multiple most parsimonious reconciliations, and chooses the reconciliation that is most consistent with neighboring families. It then re-runs makeIslands.
  
* ``printAnalysis`` produces a number of analysis/output files intended for the end user.

* ``createIslandBed`` produces bed files for each genome.

Locus families and locus islands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A brief illustration will allow us to define some terminology used in xenoGI's output. The basic goal of xenoGI is to group genes with a common origin and map them onto a phylogenetic tree.

Consider a clade of three species: (A,B),C. In this group, A and B are most closely related, and C is the outgroup. Gene a in species A has an ortholog b in species B. These two genes have high synteny, but have no ortholog in C. We call a and b a *locus family* because they are descended from a common ancestor, and occur in the same syntenic location.

When a genomic island inserts as a part of a horizontal transfer event, it typically brings in multiple locus families at the same time. xenoGI will attempt to group these into a *locus island*. In the a/b case, if there were several other locus families nearby that also inserted on the branch leading to the A,B clade, we would group them together into a single locus island.

Initial families, origin families and the DTLOR model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In fact, a locus family has several possible origins. It may be due to a horizontal transfer event coming from some other genome. Alternatively, it may reflect a rearrangement event within a genome, moving genes to a new syntenic location (for example in conjunction with a duplication event). A final possibility is that it is a core family and originated in the common ancestor of the strains under consideration. One of xenoGI's goals is to distinguish between these possibilities for each locus family (and also for the locus islands that contain them).

xenoGI does this during the process of family formation. It begins by forming large gene groupings using single linkage clustering and sequence similarity as determined by blast. It then takes these "blast families", breaks up the larger ones (which must be done for reasons of time efficiency in later steps), and uses them as a basis for making a set of families which we call initial families. For each initial family, xenoGI creates a gene tree using MUSCLE and FastTree (the user can determine whether this should be done with DNA or protein by setting the input parameter dnaBasedGeneTrees). It then reconciles each resulting gene tree to the species tree using the DTLOR model.

DTLOR is an extension we have developed to the DTL (duplication-transfer-loss) reconciliation model. It is especially suited to reconciliation in clades of closely related microbes because it allows some of the evolution of a gene family to occur outside of the given species tree. In particular, it allows multiple entry events into the species tree (where DTL allows only one). To facilitate the recognition of such entry events, the model also keeps track of the *syntenic region* of each gene as it evolves in the species tree. Two genes are said to be in the same syntenic region if they share a substantial fraction of core genes in a relatively large window around them and, second, they share a certain amount of similarity among all genes in a smaller window around them. Thus, in addition to duplication, transfer, and loss events, the DTLOR model adds *origin* events to indicate that a gene is transferred from outside of the species tree and *rearrangement* events that account for changes in the syntenic regions of genes within the same the genome.

xenoGI obtains a reconciliation for each initial family, and then uses these to break the initial families up according to origin events. The new families that result from this are called *origin families* because each one has an origin event at its base. Origin events can either correspond to core genes (if they occur at the root of the species tree) or to horizontal transfer events (if they occur below the root). In general, users will be more interested in origin families than initial families. However the class representing initial families does contain some information (the raw reconciliation output) which isn't present in the origin families, and may occasionally be of interest.

It may be helpful to give an example of the sort of thing one might find in an origin family. Consider a clade of four species: ((W,X),Y),Z::

              _____ W
         ____|s2
    ____|s1  |_____ X
   |    |
  _|s0  |__________ Y
   |
   |_______________ Z

We've labeled the internal nodes on this tree s0,s1, and s2.

Imagine that genes w1 and x1 represent a locus family in the W,X clade. They are orthologs sharing high synteny. (And they have no ortholog in species Y or Z). Imagine that there is also a paralog x2 that occurs in a different syntenic region (and that there is no w2, y2 or z2, ie W, Y and Z have no paralogs in this syntenic region). This situation could arise if there had been a horizontal transfer from outside the clade on the lineage leading to s2, and then a subsequent duplication and rearrangement after s2 on the lineage leading to X. If this were the case, xenoGI would place x1, y1, and x2 into a single origin family. w1 and x1 would be put in one locus family, and x2 in another. (In general, an origin family consists of one or more locus families.)
  
Notes on several input parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``rootFocalClade`` defines the focal clade where we will do the reconstruction. It is specified by giving the name of an internal node in the species tree. It should be chosen such that there are one or more outgroups outside the focal clade. These outgroups help us to better recognize core genes given the possibility of deletion in some lineages. 

* ``numProcesses`` determines how many separate processes to run in parts of the code that are parallel. If you have a machine with 32 processors, you would typically set this to 32 or less.

* ``dnaBasedGeneTrees`` specifies what will be used to make gene trees. If this is set to True, the method will use DNA based alignments, otherwise it will use protein alignments.

* ``useGeneRaxToMakeSpeciesTrees``. If set to True, xenoGI uses GeneRax in addition to FastTree to make species trees. GeneRax produces species-tree-aware gene trees, which are known to be of higher quality than gene trees calculated from gene sequences alone. (The cost is that GeneRax is slower). If using GeneRax then you also need to specify the parameter ``geneRaxPath``.
  
* The DTLOR cost parameters: ``duplicationCost``, ``transferCost``, ``lossCost``, ``originCost``, ``rearrangeCost``. The parsimony based reconciliation algorithm finds the minimum cost mapping of a gene tree onto the species tree. These parameters specify the costs for each of the DTLOR operations. The params.py file included in the example directory contains a set of costs we have found to work reasonably well, however users may potentially want to adjust these. The same parameters are used for all reconciliations, with one exception (see next bullet).

* ``reconcilePermissiveOriginGeneListPath``. This parameter is commented out by default, and will only be useful in certain situations. There are some genomic islands that insert repeatedly in the same syntenic region. An example is the SCCmec element in *Staphylococcus aureus*. In such cases, it is desirable to do the reconciliation with cost parameters that are permissive to origin events. xenoGI allows users to identify families that should be handled in this way. The first step is to create a file of xenoGI genes belonging to such families (one gene per line). We then set the ``reconcilePermissiveOriginGeneListPath`` to point to this file. The script ``getProteinsWithBlastHitsVsMultifasta.py`` in the misc/ directory may be useful in producing this file. The documentaiton for the misc directory has some further information.

Output files
~~~~~~~~~~~~

The last two steps, printAnalysis and createIslandBed make the output files relevant to the user.

* ``printAnalysis``

  - This script produces a set of species specific genome files. These files all have the name ``genes`` in their stem, followed by the strain name, and the extension .tsv. In the example/ data set, ``genes-E_coli_K12.tsv`` is one such. These files contain all the genes in a strain laid out in the order they occur on the contigs. Each line corresponds to one gene and contains:
    + gene name
    + origin of the gene, specified by a single character: a C indicating core gene, or an X indicating xeno horizontal transfer. This field is an interpretation of the O event from the DTLOR reconcilation based on its placement in the species tree.
    + gene history, specified by a string. This gives the history of the gene from its origin until the tip of the gene tree, and consists of single letters corresponding to the operations in the reconcilation model. D, duplication; T, transfer (within the species tree); O, origin; R, rearrangement; S, cospeciation.
    + locus island number
    + initial family number
    + origin family number
    + locus family number
    + gene description

  - ``islands.tsv`` tab delimited listing of locus islands. Each line corresponds to one locus island. The first field is the locus island number, the second field is its mrca (most recent common ancestor), and the third is a string giving the origin of each locus family in the locus island (possible values for each locus family are C for core gene, X for xeno HGT, and R for rearrangement). Subsequent fields give the locus families in this locus island. Each locus family is listed with its number, and then the genes it contains, separated by commas.
  
  - ``islandsSummary.txt`` A more human readable summary of locus islands, organized by node. This includes a tabular printout of the island, as well as a listing of each gene and its description if any.

* ``createIslandBed`` creates a subdirectory called bed/ containing bed files for each genome showing the locus islands in different colors. (Color is specified in the RGB field of the bed).

Interactive analysis
~~~~~~~~~~~~~~~~~~~~

After you have done runAll, it is possible to bring up the interpreter for interactive analysis::

  xenoGI params.py interactiveAnalysis
  
From within python, you can then run functions such as

* printLocusIslandsAtNode

  Usage::

    printLocusIslandsAtNode('s2')         # All locus islands at node s2
    printLocusIslandsAtNode('E_coli_K12') # All locus islands on the E. coli K12 branch

* findGene

  Usage::
  
    findGene('gadA')

  Find information about a gene. Searches all the fields present in the geneInfo file, so the search string can be a locus tag, protein ID, a common name, or something present in the description. For each hit, prints the gene, LocusIsland, initialFamily, originFamily, LocusFamily and gene description.
  
* printLocusIsland

  Say we've identified locus island 1550 as being of interest. We can print it like this::

    printLocusIsland(1550,10) # First argument is locus island id, second is the number of genes to print to each side
    
  printLocusIsland prints the locus island in each strain where it's present. Its output includes the locus island and family numbers for each gene, the most recent common ancestor (mrca) of the family, and a description of the gene.

* printFam

  Print scores within a particular gene family, and also with similar genes not in the family::
  
    printFam(originFamiliesO,5426)

  This function also prints a summary of the reconciliation between the gene tree for this family and the species tree.
    
  Note that this function takes a family number, not a locus family number.

Obtaining a species tree if you don't already have one
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Having an accurate species tree is a key to the xenoGI method.

The package does include some functions that may be helpful if you don't have a species tree. These use MUSCLE and FastTree to make gene trees, and ASTRAL to consolidate those gene trees into a species tree.

You begin by running the first three steps of xenoGI::

  xenoGI params.py parseGenbank
  xenoGI params.py runBlast
  xenoGI params.py calcScores

You can then run ``makeSpeciesTree``::

  xenoGI params.py makeSpeciesTree

In the ``params.py`` file, the parameter ``dnaBasedGeneTrees`` determines whether DNA or protein are used to make genes trees. (If True, DNA is used).

In order to use ``makeSpeciesTree``, you will also need to add one parameter to ``params.py``. There should be a parameter outGroup which specifies a single outgroup species to be used in rooting the species tree.

Once ``makeSpeciesTree`` has completed, you can proceed with the rest of xenoGI::

  xenoGI params.py makeFamilies
  xenoGI params.py makeIslands
  xenoGI params.py refine
  xenoGI params.py printAnalysis
  xenoGI params.py createIslandBed
  
Additional flags
~~~~~~~~~~~~~~~~

Print the version number::
   
  xenoGI params.py version

Calculate the amino acid identity between strains::

  xenoGI params.py aminoAcidIdentity

This uses blast output, and so should be run after the runBlast step. It identifies the best reciprocal hits between each pair of strains. It then averages protein identity across these, weighted by alignment length.
  
Produce a set of pdf files showing histograms of scores between all possible strains::

  xenoGI params.py plotScoreHists
  
    
Additional files
----------------

The github repository also contains an additional directory called misc/. This contains various python scripts that may be of use in conjunction with xenoGI. Installation via pip does not include this, so to use these you need to clone the github repository. There is some brief documentation included in the misc/ directory.
