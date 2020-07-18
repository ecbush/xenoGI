===============================================
Brief documentation for helper scripts in misc/
===============================================

The scripts described here are not part of the main xenoGI functionality, and can't be run with the normal "xenoGI params.py flag" syntax. To use these scripts, you need to have downloaded a copy of the repository from GitHub. You can then run the scripts by calling python on them directly as described below.


Tools for visualization with the IGB browser
--------------------------------------------

* With the IGB browser (http://bioviz.org/igb/), the bed file option allows you to most easily display different islands in different colors.

  As described in README.rst, bed files can be created via::

    xenoGI params.py createIslandBed
           
* To use with IGB , we've included scripts for making an IGB quickload directory (in misc).

  These require some additional sequence files from NCBI. For the example, those can be downloaded by going to the ncbi/ directory and running:::

  sh getSeqs.sh

  Then in the main directory, run::

    python3 path-to-xenoGI-github-repository/misc/createIgbQuckloadDirs.py ncbiIgbDirMap.txt ncbiHumanMap.txt ncbi/ igbExample

  This script requires two programs from the blat suite, faToTwoBit and twoBitInfo (https://genome.ucsc.edu/goldenpath/help/blatSpec.html). These are required to be in the path. On windows it may be easier to simply edit createIgbQuckloadDirs.py, entering the absolute path to these executables.

  Move species.txt and contents.txt into the directory igbExample.

  Now, running the script::

    sh moveInBed.sh

  moves the bed files into the newly created igbExample directory. You can then set IGB up to load this.

Obtaining alignments
--------------------

The ``makeSpeciesTree`` flag described in the main xenoGI README calculates a gene tree for each set of "all around best reciprocal hit" orthologs (ie the hard core). In the process it creates a set of either DNA or protein alignments using MUSCLE and then calls FastTree on these. If you would like to skip the FastTree step and only obtain alignments, you can use the following::

  python3 path-to-xenoGI-github-repository/misc/createAlignments.py params.py "dna" aabrhHardCore.out alignDir

If you want protein alignments rather than DNA, then change "dna" to "prot" in the above.
  
Making a tree suitable for xenoGI
---------------------------------

If you have a species tree for your strains obtained from elsewhere (ie from the phylogenetic reconstruction program of your choice) you will likely have to make several modifications to it before in can be used with xenoGI. The prepareSpeciesTree.py script in misc/ will take an input tree in Newick format. It then roots this tree, adds names to its internal nodes, and removes branch lengths.

Run it like this::

  python3 path-to-xenoGI-github-repository/misc/prepareSpeciesTree.py input.tre output.tre Vibrio_cholerae_N16961

The final argument here is a known outgroup. It can also be a clade, in which case you would separate the names of the strains in that clade with whitespace.


Reconcile a single gene tree with the species tree
--------------------------------------------------

The reconcileOneGeneTree.py script in misc/ allows you to reconcile a single gene tree against the species tree. This may be useful if you want to play with the DTLOR parameters to see how they affect the reconcilation. (To be honest it may be of most interest to developers).

This script depends on xenoGI having been run up to the makeFamilies stage, because it makes use of the syntenic location assignments that get done at that point.

As arguments it takes a species tree, a gene tree (both in newick) a xenoGI parameters file, and the costs for each of the DTLOR operations.

Run it like this::

  python3 path-to-xenoGI-github-repository/misc/reconcileOneGeneTree.py example.tre geneFamilyTrees/initFam001699.tre params.py 1 1 1 1 1

As output it produces a rooted gene tree (determined by trying all possible rootings and picking the one with the best reconcilation) and a representation of the reconcilation given by traversing the gene tree.
