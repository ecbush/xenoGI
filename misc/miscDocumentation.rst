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


Identify proteins with similarity to provided multifasta
--------------------------------------------------------

The ``getProteinsWithBlastHitsVsMultifasta.py`` script in misc/ identifies xenoGI proteins with significant similarity to a provided protein multifasta, printing their xenoGI gene names to standard out, one per line. It expects to be run in a xenoGI working directory.

Run it like this::
  
  python3 path-to-xenoGI-github-repository/misc/getProteinsWithBlastHitsVsMultifasta.py params.py strainInfo.txt fasta/multiFastaWithProtsToSearch.fa > listOfHits.txt

In a typical xenoGI working directory, ``strainInfo.txt`` is a file that should already exist listing all the strains. If you wish to run this script in an xlMode directory, and only run on strains in the scaffold, then you would first need to produce a strainInfo.txt-like file with only the scaffold strains. This could be done by running ``listTreeStrains.py`` (described below) on the scaffold tree.
  
Here's an example of how this script might be used.

As described in the README, it is possible to provide xenoGI with a list of genes on which we should use DTLOR parameters that are permissive to origin events. A case where one might want to do this would be the SCCmec element in *Staphylococcus aureus*. Say we're running xenoGI on a set of *S. aureus* genomes. We could do as follows. First, we collect a set of protein sequences from known SCCmec elements (e.g. from NCBI). Put these in a file (``fasta/multiFastaWithProtsToSearch.fa``. Then run the ``getProteinsWithBlastHitsVsMultifasta.py`` script on it from the top level of the xenoGI working directory. If we run it as described above, the output will go in to a file ``listOfHits.txt`` at the top level of the xenoGI working directory. We can now add the following line to ``params.py``::

  reconcilePermissiveOriginGeneListPath = 'listOfHits.txt'

Now when we run xenoGI DTLOR will use permissive origin costs for all families with genes in the hit list.

Get a list of the strains in a tree file
----------------------------------------

The ``listTreeStrains.py`` script in misc/ takes a parameter file and a tree as input, and produces a listing of the strains in that tree::

  python3 path-to-xenoGI-github-repository/misc/listTreeStrains.py xlParams.py scaffold.tre > scaffoldStrains.txt

xlMode: running on larger numbers of strains
--------------------------------------------

The ``runXlMode.py`` script allows you to run on many hundreds of strains. The basic strategy it follows is:

- Create a species tree for the whole data set
- Pick a subset of strains to carry out further analysis on. We call this the scaffold. The users specifies the number of strains that should be in the scaffold, and then these are chosen to maximize branch length (ie get he maximum amount of diversity). Alternatively the user can directly specify the scaffold strains.
- Run regular xenoGI on the scaffold
- Map genes for the full data set back onto the scaffold, and use this to assign them to families.

In the end this produces a xenoGI analysis on the scaffold, and a mapping so that for any gene in the whole data set, you can see what scaffold family (if any) it has been assigned to.

To run it, you would create a directory for this analysis. Inside set up an ncbi subdirectory with gbff files, just as in regular xenoGI. Copy ``xlParams.py`` from the ``misc/`` directory of the repository into this directory.

You will need to edit a few parameters in this file.

Set the ``outGroup`` parameter to specify the name of the assembly you will use as outgroup. (This makes it possible to root the species tree).

Set ``trimLeafNum`` to specify the number of strains you want in the scaffold tree. If you want to directly specify the strains to be included in the scaffold, the set ``userSpecifiedStrainsFileName`` to point to a strain file. (This strain file should contain one strain per line, and the number of strains in it must be less than trimLeafNum).

Scaffold formation involves an initial step to get a preliminary scaffold, followed by a second refinement step where some additional strains are added. The parameter ``numStrainsToAddToScaffold`` specifies how many strains to add to the scaffold for this second iteration.

To run ``runXlMode.py`` first parse the genbank files::

  python3 path-to-xenoGI-github-repository/misc/runXlMode.py xlParams.py parseGenbank

Then create sets of orthologs from core genes (to be used in tree reconstruction)::
  
  python3 path-to-xenoGI-github-repository/misc/runXlMode.py xlParams.py obtainCoreOrthoSets

Make the species tree (using MUSCLE, FastTree, and Astral)::
  
  python3 path-to-xenoGI-github-repository/misc/runXlMode.py xlParams.py makeSpeciesTree

Create the scaffold tree::
  
  python3 path-to-xenoGI-github-repository/misc/runXlMode.py xlParams.py makeScaffold

Map all genes onto the scaffold and create output files::
  
  python3 path-to-xenoGI-github-repository/misc/runXlMode.py xlParams.py printAnalysisXL

Output files can be found in the analysis/ directory. There are the normal xenoGI output files, such as genes files, ``islandsSummary.txt``, ``islands.tsv``. The file ``xlAnalysisSummary.txt`` gives the number of all genes which map onto the scaffold, and the number which does not.

If you want to know what locus family a particular gene has been assigned to (when the genes were mapped to the scaffold) then you can enter interactive more::

  python3 path-to-xenoGI-github-repository/misc/runXlMode.py xlParams.py interactiveAnalysis

Say you wanted to know the locus family number for gene 50343.
At the python prompt, you would type::

  >>> geneToLocFam(50343)
  44

If it is unmapped, this function will return None.

If you subsequently want to know more about locus family 44, you could enter interactive mode for regular xenoGI (which will work on the scaffold)::

  python3 path-to-xenoGI-github-repository/xenoGI-runner.py xlParams.py interactiveAnalysis

Then get that locus family::

  >>> lfO = originFamiliesO.getLocusFamily(44)

Get the origin family it comes from::
  >>> famNum = lfO.famNum

And print info on this origin family::

  >>> printFam(famNum,originFamiliesO)


downloadGenbank.py: automatically setting up a xenoGI working directory
-----------------------------------------------------------------------

Code for downloading ncbi genbanks. The below command is intended to
be run before any calls to xenoGI, as it sets up the ncbi folder and
human-readable map file, which are required for parseGenbank. The name
of the folder and the map file are read in from params.py, so
params.py is required, and any changes to params.py should be made
before running this script.

In addition, downloadGenbank is dependent on ncbiPythonTools.py, so that module must 
also be accessible

Calling the Function
~~~~~~~~~~~~~~~~~~~~

downloadGenbanks.py takes in three commandline arguments.
These arguments are:

1) A file name. This should be a text file with no header, and each line
   represents a genbank to be downloaded. The first column of this file
   is made of UID numbers and/or ascenscion numbers. All numbers must be from 
   the same database, which is specified as the second argument.
   The second column, separated by a tab character, is optional. This second 
   column should contain names to be associated with the genbank, which are used
   in the human map file (giving human-readable names to the downloaded genbanks)

   In the below example file, the id numbers are from the nucleotide database

   ::

      NZ_CP009044.1   
      NZ_CP007773.1   user_inputted_name
      NZ_CP020478.1   
   
   The next example file, with accession numbers, would yield the same results

   ::

      GCF_000736415.1   some_other_name
      GCF_000816305.1   
      GCF_002080395.1


2) The database associated with the id numbers. In the first file above, this would be
   ‘nuccore’, as nuccore is the name of the Entrez database. In the second file above, 
   either 'nuccore' or 'assembly' would be accepted. All possible E-utility database 
   names are listed here:
   https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly

3) an email: Entrez uses an email to make calls. While it runs
   without it, providing an email is preferred for working with the ncbi database.

Example call::
  python3 path-to-xenoGI-github-repository/misc/downloadGenbank.py assemblyList.txt assembly researcher@hmc.edu

Functions in Interactive mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If no inputs are given, downloadGenbank will start as an interactive. 
While reading in from a file is still possible using the function fileToDownload, 
other functions are included for other input types. For these functions, the email 
(if so desired) should be manually set, using ``Entrez.email = yourEmail@email.address``

* queryToDownload

  This function takes in a keyword(s) as a string, a number of results
  to be returned as an integer (retmax), and an email. Using Entrez's esearch, the specified database is
  searched using the keyword / keywords, and the search results are downloaded. 
  More on Entrez's esearch can be found here: (https://biopython.org/docs/1.75/api/Bio.Entrez.html#Bio.Entrez.esearch)

* searchToDownload

  If a more specific search is desired, an already-created Entrez search handle can be inputted.
  This function takes in the search and the database the search targets, and an email.

* downloadMultipleGBFFS

  This function is called by the other two, and takes in a python list of ids, 
  the database they are related to, and an email.

