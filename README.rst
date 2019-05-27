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

Installation
------------

The easiest way to install is using pip::

  pip3 install xenoGI


How to use
----------

An ``example/`` directory is included in this repository.

We work with a set of species with known phylogenetic relationships. In the example, these species are: E. coli K12, E. albertii, E. fergusonii and S. bongori.

Required files
~~~~~~~~~~~~~~

The working directory must contain:

* A parameter file. In the provided ``example/`` directory this is called ``params.py``. The ``blastExecutDirPath`` parameter in this file should be edited to point to the directory where the blastp and makeblastdb executables are.

* A newick format tree representing the relationships of the strains. In the example this is called ``example.tre``. Note that branch lengths are not used in xenoGI, and ``example.tre`` does not contain branch lengths. Also note that internal nodes should be given names in this tree. In the example.tre we label them i0, i1 etc. The parameter ``treeFN`` in ``params.py`` has the path to this tree file.

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
  
* ``runBlast`` does an all vs. all protein blast of the genes in these strains. The number of processes it will run in parallel is specified by the numThreads parameter in the parameter file.
  
* ``calcScores`` calculates similarity and synteny scores between genes in the strains. It is also (mostly) parallelized.
  
* ``makeFamilies`` calculates gene families in a tree aware way, also taking account of synteny.

* ``makeIslands`` groups families according to their origin, putting families with a common origin together as islands. It is partly parallelized.

* ``printAnalysis`` produces a number of analysis files.

* ``createIslandBed`` produces bed files for each genome.
  

Output files
~~~~~~~~~~~~

The last two steps, printAnalysis and createIslandBed make the output files relevant to the user.

* ``printAnalysis``

  ``islandsSummary.out`` contains a summary of islands, organized by node.

  This script also produces a set of species specific genome files. These contain all the genes in a strain laid out in the order they occur on the contigs. Each gene entry includes island and family information, as well as a brief description of the gene's function. These files all have the name genes in their stem, followed by the strain name, and the extension .out.

* ``createIslandBed`` creates a subdirectory called bed/ containing bed files for each genome showing the islands in different colors. (Color is specified in the RGB field of the bed).

Interactive analysis
~~~~~~~~~~~~~~~~~~~~

After you have done runAll, it is possible to bring up the interpreter for interactive analysis::

  xenoGI params.py interactiveAnalysis
  
From within python, you can then run functions such as

* printIslandsAtNode

  Usage::

    printIslandsAtNode('i0')         # All islands at node i0
    printIslandsAtNode('E_coli_K12') # All islands on the E. coli K12 branch

* findIsland

  Usage::
  
    findIsland('gadA') # Find an island associated with a gene name or description``
    
* printIsland

  If we've identified an island of interest (for example island number 3500) then we can print it like this::

    printIsland(3500,10) # First argument is island id, second is the number of genes to print to each side
    
  printIsland prints the island in each strain where it's present. Its output includes the island and family numbers for each gene, an error score for the family of each gene, the most recent common ancestor (mrca) of the family, and a description of the gene. The error score is intended to indicate confidence in the correctness of the family. 0 means more confident, higher numbers less confident.

* printFam

  Print scores within a particular gene family, and also with similar genes not in the family::
  
    printFam(3500)


Additional flags
~~~~~~~~~~~~~~~~

::
   
  xenoGI params.py version

will print the version number, and::


  xenoGI params.py plotScoreHists

will produce a set of pdf files showing histograms of scores between all possible strains.
  
    
Additional files
----------------

The github repository also contains an additional directory called misc/. This contains various python scripts that may be of use in conjunction with xenoGI. Installation via pip does not include this, so to use these you need to clone the github repository. There is some brief documentation included in the directory.
