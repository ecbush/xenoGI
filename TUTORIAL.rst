===============
xenoGI tutorial
===============


Introduction
------------

Before doing this tutorial, you should go through the README file, and run the example data set that is distributed with the repository. To do this you will also need to set up the various software packages that xenoGI makes use of (see the README in the repository for details).

The purpose of this tutorial is to provide an example of how one would download new data and run it in xenoGI, and also of how to look at the output. It is intended to be run on Linux or Mac. (But could probably work on Windows with some modifications).

We will assume some familiartity with the unix command line, but will otherwise try to work at a basic level.

For consistency in terms of file paths we will set up a tutorial directory and place a copy of the xenoGI repository, as well as the data set we're working on within it. To do this, start a terminal on the machine you're going to work on. In whatever location suits you, create the following directory::

  mkdir xgiTutorial

Then move into it::

  cd xgiTutorial

Next we're going to put a copy of the repository inside this directory. (Note that if you prefer to use a pip-installed version you can do that, and skip this step. Just make sure the pip-installed version corresponds to the version for this tutorial. See below). If you already have a copy, feel free to move it here (e.g. with the Mac gui, with the ``mv`` command etc). Alternatively you can clone a new copy. You can do this with a browser via the github website (https://github.com/ecbush/xenoGI). You can also do it at the command line (assuming you have git installed)::

  git clone https://github.com/ecbush/xenoGI

Next put the repository into the right release::

  cd xenoGI
  git checkout v3.0.0
  cd ..

First steps
-----------

Getting the genomes
~~~~~~~~~~~~~~~~~~~

We'll also create a directory for the data set we're going to work with::

  mkdir enterics
  cd enterics

(You should now be in xgiTutorial/enterics).
  
xenoGI makes heavy use of synteny, and therefore requires genomes that are assembled to the scaffold level or better. We're going to work with a set of complete assemblies from enteric bacteria (this is a different set than the one distributed with the repository in the example/ directory)

Within the enterics directory we'll create a subdirectory to hold the gbff files we're going to get::

  mkdir ncbi
  cd ncbi

The data set for this tutorial consists of 5 bacterial assemblies::
  
  GCF_000006945.2
  GCA_000018625.1
  GCA_000439255.1
  GCA_000005845.2
  GCF_000027085.1

If you are on linux (and have wget) you can obtain them at the command line by running this::

  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/018/625/GCA_000018625.1_ASM1862v1/GCA_000018625.1_ASM1862v1_genomic.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/439/255/GCA_000439255.1_ASM43925v1/GCA_000439255.1_ASM43925v1_genomic.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/085/GCF_000027085.1_ASM2708v1/GCF_000027085.1_ASM2708v1_genomic.gbff.gz

Alternatively, you can go to NCBI's .. _assembly https://www.ncbi.nlm.nih.gov/assembly page and enter the assembly accessions above. After clicking on the resulting link to an assembly, you should see the main page for that assembly. On the right there will be a section entitled ``Download Assembly``. You will now click on a link to the ftp directory. For GCF_000006945.2 and GCF_000027085.1 click on the link to the RefSeq assembly. For the others the GenBank assembly. On the ftp page, download the file ending gbff.gz. Once you have all of these, move them into ``xgiTutorial/enterics/ncbi``.

We next need to uncompress these assembly files, which can be done with::

  gunzip *.gz

Finally, move up one level to the enterics directory::

  cd ..

(You should now be in ``xgiTutorial/enterics/``).
  
Setting up a few required files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy over a parameter file from the example directory::

  cp ../xenoGI/example/params.py .

Next edit ``params.py`` using a text editor such as emacs, vim, nano etc. You should edit the following to give the absolute (full) path to the directory where the blastp and makeblastdb executables reside::

  blastExecutDirPath = '/usr/bin/'

Also make sure that ``astralPath``, ``musclePath``, ``fastTreePath``, and ``javaPath`` all correctly point to the executables for ASTRAL, MUSCLE, FastTree and Java respectively.

Go to the ``speciesTreeFN`` entry, and edit it to read::

  speciesTreeFN='enterics.tre'

Here we've just changed the expected name of the tree file. Of course it doesn't really matter what we name the tree, but it's nice for it to correspond to the species we are working with. (Note that this tree file does not yet exist. We'll create it below).

We also need to specify the names we'll use to refer to each species. Using a text editor, create the file ``ncbiHumanMap.txt`` and paste the following into it::

  GCF_000006945.2_ASM694v2_genomic.gbff	S_enterica_LT2
  GCA_000018625.1_ASM1862v1_genomic.gbff	S_enterica_AZ
  GCA_000439255.1_ASM43925v1_genomic.gbff	S_bongori
  GCA_000005845.2_ASM584v2_genomic.gbff	E_coli_K12
  GCF_000027085.1_ASM2708v1_genomic.gbff	C_rodentium


Using screen
~~~~~~~~~~~~

xenoGI is a command line program that sometimes can take a while to run. If you are working on a remote machine, it may be useful to run xenoGI from within ``screen``, which is available on most linux distributions. For this tutorial, ``screen`` shouldn't be necessary because everything runs in a few minutes. But if you move on to larger datasets it might be helpful.

What screen does is provide a command line which you can "detach". You can then logout of the machine, and your process will keep running. When you log back in, you can retrieve it.

To start screen the first time::

  screen

To detach once you have something running::

  Ctrl+a d

To retrieve the previous screen session::

  screen -r

And finally, when you are all done and want screen to go away::

  Ctrl+d

Running xenoGI
--------------
 
Parsing files, blast etc.
~~~~~~~~~~~~~~~~~~~~~~~~~

Creating families and 

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this step we'll determine the species tree for the strains we're looking at. When working on your own data, if you already know the tree, then you would typically skip this step.

(If you don't wan't to do this step in the tutorial, you can skip to the end of this section where the correct species tree is printed, and proceed from there.)

This is the first step where various other software packages are used. It uses MUSCLE and FastTree to make gene trees, and ASTRAL to consolidate those gene trees into a species tree.

We do require that the user specify an outgroup so that we can root the species tree. In the enteric data set we're using, C_rodentium is the outgroup. Before we run the step, we need to specify the outgroup in the ``params.py`` file. Open that file in a text editor. In the 'Making species trees' section there is a parameter ``outGroup`` which has been commented out. Uncomment this (delete the hash) and set it so it reads::

  outGroup = 'C_rodentium'

Then run like so::

  python3 ../xenoGI/xenoGI-runner.py params.py makeSpeciesTree

This should take a minute or so. It will produce a newick file called ``enterics.tre``. If you skipped this step, you should manually create an ``enterics.tre`` file, with the following contents::

  ((E_coli_K12,(S_bongori,(S_enterica_LT2,S_enterica_AZ)s3)s2)s1,C_rodentium)s0;

For your reference, here's an ascii drawing of the tree, with internal nodes labelled::

         _____ E_coli_K12
    ____|
   |    |s1    ____ S_bongori
   |    |_____|
  _|          |s2   _____ S_enterica_LT2
   |s0        |____|s3
   |               |_____ S_enterica_AZ
   |
   |____ C_rodentium

  
Creating gene families and locus islands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Creating output files
~~~~~~~~~~~~~~~~~~~~~

We can now create a set of output files which we'll use in subsequent analysis::

  python3 ../xenoGI/xenoGI-runner.py params.py printAnalysis

This step is very quick, taking just a few seconds on this data set.

Examining those output files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The above command creates a subdirectory called analysis. Let's have a look at it::

  cd analysis/
  ls

You should see a set of files beginning with "genes", as well as ``islandsSummary.txt`` and ``islands.tsv``.

The genes files contain all the genes in a strain laid out in the order they occur on the contigs (the first line of each specifies what the columns are). Let's start out by looking at a known pathogenicity island, Salmonella Pathogenicity Island 1 (SPI1). This island is known to be present in all three Salmonella strains, S_enterica_LT2, S_enterica_AZ, and S_bongori. In S_enterica_LT2 it is known to extend from STM2865 to STM2900. Let's take a look::

  less -S genes-S_enterica_LT2.tsv

The -S tells the text viewer less not to wrap lines, which makes it a little easier to read. You may want to maximize your window, or make it wider so that more of each line displays. At the right of each line is included a description of each gene.

You can now search within less by entering the locus tag for the first gene STM2865.

Here's a truncated bit of what you should see::
  
  21087_S_enterica_LT2-STM2863    C       OSSS    3229    2724    3173    3229    s0      sitC - iron ABC transporter
  21088_S_enterica_LT2-STM2864    C       OSSS    3228    2723    3172    3228    s0      sitD - iron ABC transporter
  21089_S_enterica_LT2-STM2865    X       OS      3646    4170    5007    5073    s2      avrA - putative inner membr
  21090_S_enterica_LT2-STM2866    X       OSS     3646    3228    3799    3861    s2      sprB - transcriptional regu
  21091_S_enterica_LT2-STM2867    X       OSS     3646    3053    3588    3649    s2      hilC - AraC family transcri
  21092_S_enterica_LT2-STM2868    X       OSS     3646    3317    3912    3976    s2      type III secretion system e
  21093_S_enterica_LT2-STM2869    X       OSS     3646    3316    3911    3975    s2      orgA - invasion protein Org
  21094_S_enterica_LT2-STM2870    X       OSS     3646    3315    3910    3974    s2      putative inner membrane pro
  21095_S_enterica_LT2-STM2871    X       OSS     3646    3209    3770    3832    s2      prgK - EscJ/YscJ/HrcJ famil
  21096_S_enterica_LT2-STM2872    X       OSS     3646    3314    3909    3973    s2      prgJ - type III secretion s
  21097_S_enterica_LT2-STM2873    X       OSS     3646    3313    3908    3972    s2      prgI - EscF/YscF/HrpA famil
  21098_S_enterica_LT2-STM2874    X       OSS     3646    3312    3907    3971    s2      prgH - type III secretion s
  21099_S_enterica_LT2-STM2875    X       OSS     3646    3051    3586    3646    s2      hilD - AraC family transcri
  21100_S_enterica_LT2-STM2876    X       OSS     3646    3248    3827    3889    s2      hilA - transcriptional regu
  21101_S_enterica_LT2-STM2877    X       OSS     3646    3188    3748    3810    s2      iagB - invasion protein Iag
  21102_S_enterica_LT2-STM2878    X       OSS     3646    3311    3906    3970    s2      sptP - pathogenicity island
  21103_S_enterica_LT2-STM2879    X       OSS     3646    3310    3905    3969    s2      sicP - chaperone protein Si
  21104_S_enterica_LT2-STM2880    X       OS      4782    3941    4716    4782    s3      putative cytoplasmic protei
  21105_S_enterica_LT2-STM2881    X       OSS     3646    3160    3712    3774    s2      iacP - putative acyl carrie
  21106_S_enterica_LT2-STM2882    X       OSS     3646    3309    3904    3968    s2      sipA - pathogenicity island
  21107_S_enterica_LT2-STM2883    X       OSS     3646    3308    3903    3967    s2      sipD - cell invasion protei
  21108_S_enterica_LT2-STM2884    X       OSS     3646    3307    3902    3966    s2      sipC - pathogenicity island
  21109_S_enterica_LT2-STM2885    X       OSS     3646    3306    3901    3965    s2      sipB - pathogenicity island
  21110_S_enterica_LT2-STM2886    X       OSS     3646    3187    3747    3809    s2      sicA - CesD/SycD/LcrH famil
  21111_S_enterica_LT2-STM2887    X       OSS     3646    3126    3668    3730    s2      spaS - EscU/YscU/HrcU famil
  21112_S_enterica_LT2-STM2888    X       OSS     3646    3208    3769    3831    s2      spaR - EscT/YscT/HrcT famil
  21113_S_enterica_LT2-STM2889    X       OSS     3646    3207    3768    3830    s2      spaQ - EscS/YscS/HrcS famil
  21114_S_enterica_LT2-STM2890    X       OSS     3646    3125    3667    3729    s2      spaP - EscR/YscR/HrcR famil
  21115_S_enterica_LT2-STM2891    X       OSS     3646    3305    3900    3964    s2      spaO - type III secretion s
  21116_S_enterica_LT2-STM2892    X       OSS     3646    3304    3899    3963    s2      invJ - antigen presentation
  21117_S_enterica_LT2-STM2893    X       OSS     3646    3303    3898    3962    s2      invI - type III secretion s
  21118_S_enterica_LT2-STM2894    X       OSS     3646    3058    3593    3655    s2      invC - EscN/YscN/HrcN famil
  21119_S_enterica_LT2-STM2895    X       OSS     3646    3302    3897    3961    s2      invB - type III secretion s
  21120_S_enterica_LT2-STM2896    X       OSS     3646    3124    3666    3728    s2      invA - EscV/YscV/HrcV famil
  21121_S_enterica_LT2-STM2897    X       OSS     3646    3301    3896    3960    s2      invE - SepL/TyeA/HrpJ famil
  21122_S_enterica_LT2-STM2898    X       OSS     3646    3206    3767    3829    s2      invG - EscC/YscC/HrcC famil
  21123_S_enterica_LT2-STM2899    X       OSS     3646    3300    3895    3959    s2      invF - invasion protein
  21124_S_enterica_LT2-STM2900    X       OSS     3646    3299    3894    3958    s2      invH - invasion lipoprotei
  21125_S_enterica_LT2-STM2901    X       O       4591    3864    4614    4680    S_enterica_LT2  hypothetical protei
  21126_S_enterica_LT2-STM2902    X       O       4591    3802    4525    4591    S_enterica_LT2  putative cytoplasmi

The first column consists of genes listed by their xenoGI name (the locus tag is the last part of this). xenoGI has identified a locus island that corresponds to SPI1. The number for this locus island is given in column 4, and is 3646 here. (It is possible that the numbering will be different on other machines). This locus island extends from 21089_S_enterica_LT2-STM2865 to 21124_S_enterica_LT2-STM2900 as expected. Note that in the display above, we've included a few genes on either end of the locus island.

As discussed in the README, a locus island represents a set of gene families with a common origin. In this case, it corresponds to a genomic island which is inferred to have inserted on the branch leading to s2 (the branch inserted on is given in the 8th column).

Every gene in a particular clade is either a core gene, or arose by xeno horizontal transfer (horizontal transfer from outside the clade). One of the goals of xenoGI is to determine this origin for each gene. The second column in the genes file contains this information. C stands for core, and X for xeno horizontal transfer. You can note that for SPI1, all the genes are marked X.

The third column contains a gene history string. Taking the gene 21124_S_enterica_LT2-STM290 for example (invH) the string is OSS. This reflects the history of the gene after insertion, as reconstructed by the DTLOR reconciliation. O stands for origin (in this case the xeno hgt event). And S stands for co-speciation--what happens when a speciation event occurs and both descendent lineages inherit a gene. invH is inferred to have inserted on branch s2. It then underwent co-speciation events at node s2 and node s3.

As we noted, the 4th column gives the locus island. The 5th gives the initial family number, the 6th the origin family number, and the 7th the locus family number. We'll use some of these in the examples below.
