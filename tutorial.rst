===============
xenoGI tutorial
===============


Introduction
------------

Before doing this tutorial, you should go through the README file, and also run the example data set that is distributed with the repository. You will also need to set up the various software packages that xenoGI makes use of (see the README in the repository for details).

The purpose of this tutorial is to provide an example of how one would download new data and run it in xenoGI, and also of how to look at the output. It is intended to be run on Linux or Mac. (But could probably work on Windows with some modifications).

We will assume some familiartity with the unix command line, but will otherwise try to work at a basic level.

For consistency in terms of file paths we will set up a tutorial directory and place a copy of the xenoGI repository, as well as the data set we're working on within it. To do this, start a terminal on the machine you're going to work on. In whatever location suits you, create the following directory::

  mkdir xgiTutorial

Then move into it::

  cd xgiTutorial

Next we're going to put a copy of the repository inside this directory. (Note that if you prefer to use a pip-installed version you can do that, and skip this step. Just make sure the pip-installed version is the right one). If you already have a copy, feel free to move it here (e.g. with the Mac gui, with the ``mv`` command etc). Alternatively you can clone a new copy. If you're on a Mac you could do this with a browser via the github website (https://github.com/ecbush/xenoGI). You can also do it at the command line (assuming you have git installed)::

  git clone https://github.com/ecbush/xenoGI

Then put the repository into the right release::

  cd xenoGI
  git checkout v3.0.0
  cd ..


  
Getting the genomes
-------------------

We'll also create a directory for the data set we're going to work with::

  mkdir enterics
  cd enterics

xenoGI makes heavy use of synteny, and therefore requires genomes that are assembled to the scaffold level or better. We're going to work with a set of complete assemblies from enteric bacteria (this is a different set than the one distributed with the repository in the example/ directory)

Within the enterics directory we'll create a subdirectory to hold the gbff files we're going to get::

  mkdir ncbi

Setting up a few required files
-------------------------------

Copy over a parameter file from the example directory::

  cp ../xenoGI/example/params.py .

Next edit ``params.py`` using a text editor such as emacs, vim, nano etc. You should edit the following to give the absolute (full) path to the directory where the blastp and makeblastdb executables reside::

  blastExecutDirPath = '/usr/bin/'

Also make sure that ``astralPath``, ``musclePath``, ``fastTreePath``, and ``javaPath`` all correctly point to the executables for ASTRAL, MUSCLE, FastTree and Java respectively.

Go to the ``speciesTreeFN`` entry, and edit it to read::

  speciesTreeFN='enterics.tre'

Here we've just changed the expected name of the tree file. Of course it doesn't really matter what we name the tree, but it's nice for it to correspond to the species we are working with.

ncbiHumanMap.txt

Using screen
------------

Useful if you are logging in to a remote machine.

For this tutorial, everthing goes pretty fast, so may not be necessary. But if you use later with a bigger data set.


First steps
-----------

The very first thing we'll do is have xenoGI run through these genbank files, and extract the protein annotations that we'll be using::

  python3 ../xenoGI/xenoGI-runner.py params.py parseGenbank

This should take 10-15 seconds.

Note that we wrote ``python3`` above, but on some systems you may want to write simply ``python``. Just be sure that this is calling the correct version of python, with the various necessary python packages. If you are using a pip-installed verion of xenoGI, then your command would look like this::

  xenoGI params.py parseGenbank

(You can make the equivalent adjustment for the commands to follow).

Next we do an all vs. all protein blast::

  python3 ../xenoGI/xenoGI-runner.py params.py runBlast

This will take several minutes. For the steps below we will also try to give you a sense how long it should take on the tutorial data set. Note that speed may vary somewhat on your setup, but these numbers should give you a rough idea. If you subsequently do this on a larger data set of your own, of course it will take longer.

And then we calculate various types of scores::

  python3 ../xenoGI/xenoGI-runner.py params.py calcScores

This should take about 30 seconds.
  
Determining the species tree
----------------------------

In this step we'll determine the species tree for the strains we're looking at. When working on your own data, if you already know the tree, then you would typically skip this step.

(If you don't wan't to do this step in the tutorial, you can skip to the end of this section where the correct species tree is printed, and proceed from there.)

This is the first step where various other software packages are used. It uses MUSCLE and FastTree to make gene trees, and ASTRAL to consolidate those gene trees into a species tree.

We do require that the user specify an outgroup so that we can root the species tree. In the enteric data set we're using, C_rodentium is the outgroup. Before we run the step, we need to specify the outgroup in the ``params.py`` file. Open that file in a text editor. In the 'Making species trees' section there is a parameter ``outGroup`` which has been commented out. Uncomment this (delete the hash) and set it so it reads::

  outGroup = 'C_rodentium'

Then run like so::

  python3 ../xenoGI/xenoGI-runner.py params.py makeSpeciesTree

This should take a minute or so. It will produce a newick file called ``enterics.tre``. If you skipped this step, you should manually create an ``enterics.tre`` file, with the following contents::

  ((E_coli_K12,(S_bongori,(S_enterica_LT2,S_enterica_AZ)s3)s2)s1,C_rodentium)s0;

Creating gene families and locus islands
----------------------------------------

xenoGI does its most detailed reconstruction within a focal clade, leaving one or more species as outgroups. Such outgroups help us to better recognize core genes given the possibility of deletion in some lineages. One parameter we must set is the root of the focal clade. Once again, edit the ``params.py`` file. The line defining the ``rootFocalClade`` should be as follows::

  rootFocalClade = 's2'

If it doesn't already say that, change it. This says that the focal clade will be defined by the internal node ``s2``, and corresponds to the Salmonella genus. ``C_rodentium`` and ``E_coli_K12`` will be outgroups.

We will now do a series of steps to make gene families and locus islands.

Create gene families::

  python3 ../xenoGI/xenoGI-runner.py params.py makeFamilies

This will take several minutes.
  
Next create locus islands::
  
  python3 ../xenoGI/xenoGI-runner.py params.py makeIslands

This will likely take 1-2 minutes.

Then, refine families and remake islands::

  python3 ../xenoGI/xenoGI-runner.py params.py refine

This will also take 1-2 minutes. In the refinement step, xenoGI goes back and looks at cases where there are multiple most-parsimonious reconciliations. In the previous ``makeFamilies`` step, one of these was chosen arbitrarily. Now xenoGI considers all of the possibilities, and determines which of these is optimal by examining nearby gene families. (On the logic that since these will often have a common origin, it makes sense to chose the most-parsimonious reconciliation the corresponds best to them.)


Analysis
--------

We can now create a set of output files which we'll use in subsequent analysis::

  python3.7 ../xenoGI/xenoGI-runner.py params.py printAnalysis

This step is very quick, taking just a few seconds on this data set.
