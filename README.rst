======
xenoGI
======

Code for detecting genomic island insertions in clades of microbes.

Requirements
------------

* NCBI blast+

  We need blastp and makeblastdb executables (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

* Python 3

* Package dependencies

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

  If you make use of the ``makeSpeciesTree`` flag, you will also need the following

  - MUSCLE (https://www.drive5.com/muscle/).

  - FastTree (http://www.microbesonline.org/fasttree/).

  - ASTRAL (https://github.com/smirarab/ASTRAL/).

Installation
------------

The easiest way to install is using pip::

  pip3 install xenoGI

Citation
--------

If you use xenoGI in a publication, please cite the following:

Bush EC, Clark AE, DeRanek CA, Eng A, Forman J, Heath K, Lee AB, Stoebel DM, Wang Z, Wilber M, Wu H. xenoGI: reconstructing the history of genomic island insertions in clades of closely related bacteria. BMC Bioinformatics. 19(32). 2018.

How to use
----------

An ``example/`` directory is included in this repository.

The basic method works on a set of species with known phylogenetic relationships. In the example, these species are: E. coli K12, E. coli ATCC 11775, E. fergusonii and S. bongori.

Required files
~~~~~~~~~~~~~~

The working directory must contain:

* A parameter file. In the provided ``example/`` directory this is called ``params.py``. The ``blastExecutDirPath`` parameter in this file should be edited to point to the directory where the blastp and makeblastdb executables are.

* A newick format tree representing the relationships of the strains. In the example this is called ``example.tre``. Note that branch lengths are not used in xenoGI, and ``example.tre`` does not contain branch lengths. Also note that internal nodes should be given names in this tree. In the example.tre we label them i0, i1 etc. The parameter ``treeFN`` in ``params.py`` has the path to this tree file. If a strain tree is not available, xenoGI has some accessory methods, described below, to help obtain one.

* A subdirectory of sequence files. In the example, this is called ``ncbi/``. Contained in this subdirectory will be genbank (gbff) files for the species. The parameter ``genbankFilePath`` in ``params.py`` has the path to these files.

Naming of genbank files
~~~~~~~~~~~~~~~~~~~~~~~

The system needs a way to connect the sequence files to the names used in the tree.

In the example, the sequence files have names corresponding to their assembly accession number from ncbi. We connect these to the human readable names in example.tre using a mapping given in the file ``ncbiHumanMap.txt``. This file has two columns, the first giving the name of the genbank file, and the second giving the name for the species used in the tree file. Note that the species name should not contain any dashes ("-"). In ``params.py`` the parameter ``fileNameMapFN`` is set to point to this file.

Another approach is to change the names of the sequence files to match what's in the tree. If you do this, then you should set ``fileNameMapFN = None`` in ``params.py``. (This is not necessary in the example, which is already set to run the other way).

Running the code
~~~~~~~~~~~~~~~~

If you install via pip, then you should have an executable script in your path called xenoGI.

You run the code from within the working directory. To run the example, you would cd into the ``example/`` directory. You will need to ensure that the ``params.py`` parameters file contains the  correct path to the directory with the blastp and makeblastdb executables in it. Then, the various steps of xenoGI can be run all at once like this::

  xenoGI params.py runAll

They can also be run individually::

  xenoGI params.py parseGenbank
  xenoGI params.py runBlast
  xenoGI params.py calcScores
  xenoGI params.py makeFamilies
  xenoGI params.py makeIslands
  xenoGI params.py printAnalysis
  xenoGI params.py createIslandBed

If for some reason you don't want to install via pip, then you can download the repository and run the code like this::

  python3 path-to-xenoGI-github-repository/xenoGI-runner.py params.py runAll

(In this case you will have to make sure all the package dependencies are satisfied.)

What the steps do
~~~~~~~~~~~~~~~~~

* ``parseGenbank`` runs through the genbank files and produces input files that are used by subsequent code.
  
* ``runBlast`` does an all vs. all protein blast of the genes in these strains. The number of processes it will run in parallel is specified by the numThreads parameter in the parameter file. Before running a particular comparison, runBlast checks to see if the output file for that comparison already exists (e.g. from a previous run). If so it skips the comparison.
  
* ``calcScores`` calculates similarity and synteny scores between genes in the strains. It is also (mostly) parallelized.
  
* ``makeFamilies`` calculates gene families in a tree aware way, also taking account of synteny.

* ``makeIslands`` groups families according to their origin, putting families with a common origin together as islands. It is partly parallelized.

* ``printAnalysis`` produces a number of analysis files.

* ``createIslandBed`` produces bed files for each genome.
  
Notes on several parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``rootFocalClade`` defines the focal clade where we will do the reconstruction. It is specified by giving the name of an internal node. It should be chosen such that there are one or more outgroups outside the focal clade. These outgroups help us to better recognize core genes given the possibility of deletion in some lineages. 

* ``numProcesses`` determines how many separate processes to run in parts of the code that are parallel. If you have a machine with 32 processors, you would typically set this to 32 or less.


A note on the output
~~~~~~~~~~~~~~~~~~~~

A brief illustration will allow us to define some terminology used in xenoGI's output. The basic goal of xenoGI is to group genes with a common origin and map them onto a phylogenetic tree.

Consider a clade of three species: (A,B),C. In this group, A and B are most closely related, and C is the outgroup. Gene a in species A has an ortholog b in species B. These two genes have high synteny, but have no ortholog in C. We call a and b a *locus family* because they are descended from a common ancestor, and occur in the same syntenic location.

When a genomic island inserts as a part of a horizontal transfer event, it typically brings in multiple locus families at the same time. xenoGI will attempt to group these into a *locus island*. In the a/b case, if there were several other locus families nearby that also inserted on the branch leading to the A,B clade, we would group them together into a single locus island.

At present, locus islands and locus families are the basic units of output. If you are interested in finding genomic islands that inserted on a particular branch in your tree, you would be looking for locus islands identified on that branch.

Let us define one last bit of terminology. Consider another clade of three species: (X,Y),Z. Genes x1 and y1 represent a locus family in the X,Y clade. They are orthologs sharing high synteny. (And they have no ortholog species Z). Imagine that there is also a set of paralogs x2 and y2 which resulted from a gene duplication in the lineage leading to the X,Y clade. These occur in a different syntenic location. In this case, x2 and y2 constitute another locus family. Because these two locus families descended from a common ancestor gene within the species tree, we place them in the same *family*. In general, a family consists of one or more locus families.

Output files
~~~~~~~~~~~~

The last two steps, printAnalysis and createIslandBed make the output files relevant to the user.

* ``printAnalysis``

  ``islandsSummary.out`` contains a summary of islands, organized by node.

  This script also produces a set of species specific genome files. These contain all the genes in a strain laid out in the order they occur on the contigs. Each gene entry includes locus island and family information, as well as a brief description of the gene's function. These files all have the name genes in their stem, followed by the strain name, and the extension .out.

* ``createIslandBed`` creates a subdirectory called bed/ containing bed files for each genome showing the locus islands in different colors. (Color is specified in the RGB field of the bed).

Interactive analysis
~~~~~~~~~~~~~~~~~~~~

After you have done runAll, it is possible to bring up the interpreter for interactive analysis::

  xenoGI params.py interactiveAnalysis
  
From within python, you can then run functions such as

* printLocusIslandsAtNode

  Usage::

    printLocusIslandsAtNode('i2')         # All locus islands at node i2
    printLocusIslandsAtNode('E_coli_K12') # All locus islands on the E. coli K12 branch

* findLocusIsland

  Usage::
  
    findLocusIsland('gadA') # Find a locus island associated with a gene name or description``
    
* printLocusIsland

  Say we've identified locus island 1550 as being of interest. We can print it like this::

    printLocusIsland(1550,10) # First argument is locus island id, second is the number of genes to print to each side
    
  printLocusIsland prints the locus island in each strain where it's present. Its output includes the locus island and family numbers for each gene, the most recent common ancestor (mrca) of the family, and a description of the gene.

* printFam

  Print scores within a particular gene family, and also with similar genes not in the family::
  
    printFam(3279)

  Note that this function takes a family number, not a locus family number.

Obtaining a species tree if you don't already have one
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Having an accurate species tree is a key to the xenoGI method.

The package does include some functions that may be helpful if you don't have a species tree. These use MUSCLE to create either protein or DNA alignments, FastTree to make gene trees, and ASTRAL to consolidate the gene trees into a species tree.

You begin by running the first three steps of xenoGI::

  xenoGI params.py parseGenbank
  xenoGI params.py runBlast
  xenoGI params.py calcScores

You can then run ``makeSpeciesTree``::

  xenoGI params.py makeSpeciesTree

The ``params.py`` file found in the example directory contains a number of parameters related to ``makeSpeciesTree``. Among these is ``dnaBasedSpeciesTree``. If this is True, the method will use DNA based alignments, otherwise it will use protein alignments. Once ``makeSpeciesTree`` has completed, you can proceed with the rest of xenoGI::

  xenoGI params.py makeFamilies
  xenoGI params.py makeIslands
  xenoGI params.py printAnalysis
  xenoGI params.py createIslandBed
  
Additional flags
~~~~~~~~~~~~~~~~

Print the version number::
   
  xenoGI params.py version

Produce a directory containing a gene tree for every family::

  xenoGI params.py makeGeneFamilyTrees

This uses the same methods as the makeSpeciesTree flag (but doesn't call ASTRAL).
  
Produce a set of pdf files showing histograms of scores between all possible strains::

  xenoGI params.py plotScoreHists
  
    
Additional files
----------------

The github repository also contains an additional directory called misc/. This contains various python scripts that may be of use in conjunction with xenoGI. Installation via pip does not include this, so to use these you need to clone the github repository. There is some brief documentation included in the directory.
