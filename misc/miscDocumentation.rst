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

* We also include a script for creating gff files.

  It can be run like this::
    
    python3 path-to-xenoGI-github-repository/misc/createIslandGffs.py params.py


Obtaining a tree if you don't already have one
-----------------------------------------------

If you don't have a species tree for your set of species, there is some code in misc/ to help you get alignments, which can then be used (with your phylogenetic package of choice) to get a tree.

The first step is to get the genbank files for your species, and begin the process of running xenoGI, as described in README.rst. However you only run the first three steps::

  xenoGI params.py parseGenbank
  xenoGI params.py runBlast
  xenoGI params.py calcScores

Obtaining alignments
~~~~~~~~~~~~~~~~~~~~

In the process of calculating the scores, xenoGI has created sets of "all around best reciprocal hit" orthologs. These sets have exactly one copy in each strain. We now will take the proteins for each ortholog set, and construct alignments between them. This step assumes you have muscle (https://www.drive5.com/muscle/) in your path::

  python3 path-to-xenoGI-github-repository/misc/aabrhProtAlign.py params.py alignedProts.fa

One possibility is you might want to make a tree based on the protein. In that case, you can concatenate the proteins like this::

  python3 path-to-xenoGI-github-repository/misc/concatenateAlignment.py alignedProts.fa alignedProtsConcat.fa

Alternatively, you may want to work with aligned nucleotide sequences. In that case, we download some corresponding nucleotide sequences to the gbff files in example/. From within the example/ncbi/ directory, run::

  sh getSeqs.sh
  gunzip *.gz

Now back in example/ run::
  
  python3 path-to-xenoGI-github-repository/misc/aabrhBackAlign.py 4 alignedProts.fa alignedNucs.fa ncbi/*_cds_from_genomic.fna

  python3 path-to-xenoGI-github-repository/misc/concatenateAlignment.py alignedNucs.fa alignedNucsConcat.fa

  Note that one of the inputs to aabrhBackAlign.py is the number of strains (4 in this case).
  
Making a tree suitable for xenoGI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now you have output from some phylogenetic reconstruction program in the form of a newick tree. You will very likely have to make several modifications to it before in can be used with xenoGI. misc/ contains scripts for this.

* Make it a rooted tree

  Assuming you have a file from a phylogenetic reconstruction program called unrooted.tree, run::
    
    python3 path-to-xenoGI-github-repository/misc/rootTree.py unrooted.tre rooted.tre Vibrio_cholerae_N16961
  
  The final argument here is a known outgroup. It can also be a clade, in which you would separate the names of the strains in that clade with whitespace.

* Add names to the internal nodes

  xenoGI requires internal nodes to be named. We can do that like this::

    python3 path-to-xenoGI-github-repository/misc/nameInternalNodes.py rooted.tre namedNode.tre

* Remove branch lengths (optional)::
    
    python3 path-to-xenoGI-github-repository/misc/stripBranchLen.py namedNode.tre final.tre

Once you have a tree, you can resume with the normal xenoGI steps (makeFamilies, makeIslands, printAnalysis) as described in README.rst.
