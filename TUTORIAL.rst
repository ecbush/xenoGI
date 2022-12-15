===============
xenoGI tutorial
===============


Introduction
------------

Before doing this tutorial, you should go through the README file, and run the example data set that is distributed with the repository. The README also has notes on installation and the various software packages that are needed.

The purpose of this tutorial is to provide an example of how one would download new data and run it in xenoGI, and also of how to look at the output. It is intended to be run on Linux or Mac. (But could probably work on Windows with some modifications).

We will assume some familiartity with the unix command line, but will otherwise try to work at a basic level.

First steps
-----------

Setting up the tutorial directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We'll first set up a directory for the tutorial. Call it ``xgiTutorial``. (You can do this with the command line, or using the mac gui if you're using a mac).

Next we're going to put a copy of the repository inside this directory. (Note that if you prefer to use a pip-installed version you can do that, and skip this step. Just make sure the pip-installed version corresponds to the version for this tutorial. See below). If you already have a copy of the repository, feel free to move it here. Alternatively you can clone a new copy, either with a browser via the github website (https://github.com/ecbush/xenoGI) or via the command line (assuming you have git installed)::

  git clone https://github.com/ecbush/xenoGI

Next put the repository into the right release. At the command line, cd into the repository directory (``xgiTutorial/xenoGI``). Then do::

  git checkout v3.1.1

Next create a second subdirectory of ``xgiTutorial/`` called ``enterics/``. This is the working directory for the data we'll be using.

Within ``enterics/`` create a subdirectory ``ncbi/``. This is where we're going to put the genome assemblies. (To be clear the path for this directory should be ``xgiTutorial/enterics/ncbi/``).

Getting the genomes
~~~~~~~~~~~~~~~~~~~

We'll now download the genome sequences and put them in the ``ncbi`` folder.

xenoGI makes heavy use of synteny, and therefore requires genomes that are assembled to the scaffold level or better. We're going to work with a set of complete assemblies from enteric bacteria (this is a different set than the one distributed with the repository in the example/ directory)

The data set for this tutorial consists of 5 bacterial assemblies with the following GenBank accession numbers::
  
  GCF_000006945.2
  GCA_000018625.1
  GCA_000439255.1
  GCA_000005845.2
  GCF_000027085.1

If you are on linux (and have wget) you can obtain the files at the command line by running this::

  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/018/625/GCA_000018625.1_ASM1862v1/GCA_000018625.1_ASM1862v1_genomic.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/439/255/GCA_000439255.1_ASM43925v1/GCA_000439255.1_ASM43925v1_genomic.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.gbff.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/085/GCF_000027085.1_ASM2708v1/GCF_000027085.1_ASM2708v1_genomic.gbff.gz

Alternatively, you can go to NCBI's assembly page (https://www.ncbi.nlm.nih.gov/assembly) and enter the assembly accessions above. xenoGI uses genbank gbff files as input. You'll be downloading files ending in ``gbff.gz``. For more information on how to download genome sequences from NCBI please visit https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/.

After you've downloaded the assemblies and put them in the ``ncbi/`` directory, you need to uncompress them. This can be done at the command line with::

  gunzip *.gz
  
Setting up ``params.py``
~~~~~~~~~~~~~~~~~~~~~~~~

Move back into the main working directory ``xgiTutorial/enterics/``.

We'll next obtain a parameter file by copying the one in the ``example/`` directory distributed with the repository::

  cp ../xenoGI/example/params.py .

Next edit ``params.py`` so that ``blastExecutDirPath`` correctly indicates the directory where the blastp and makeblastdb executables reside.

Also make sure that ``astralPath``, ``musclePath``,``geneRaxPath``, ``fastTreePath``, and ``javaPath`` contain the correct paths to the executables for ASTRAL, MUSCLE, GeneRaxx, FastTree and Java respectively.

Then go to the ``speciesTreeFN`` entry, and edit it to read::

  speciesTreeFN='enterics.tre'

Here we've just changed the expected name of the tree file. Of course it doesn't really matter what we name the tree, but it's nice for it to correspond to the species we are working with. (Note that this tree file does not yet exist. We'll create it below).

You may also want to edit the ``numProcesses`` entry. This determines how many separate processes to run in parts of the code that are parallel. If you have a machine with 32 processors, you would typically set this to 32 or less.

Setting up ``ncbiHumanMap.txt``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We also need to specify the names we'll use to refer to each strain. Using a text editor, create the file ``ncbiHumanMap.txt`` and paste the following into it::

  GCF_000006945.2_ASM694v2_genomic.gbff	S_enterica_LT2
  GCA_000018625.1_ASM1862v1_genomic.gbff	S_enterica_AZ
  GCA_000439255.1_ASM43925v1_genomic.gbff	S_bongori
  GCA_000005845.2_ASM584v2_genomic.gbff	E_coli_K12
  GCF_000027085.1_ASM2708v1_genomic.gbff	C_rodentium

Note that it is possible to set up human readable names by renaming the assembly files. We won't do that here, but if you want to do that for your own projects, see the README. Also see the README for restrictions on characters allowed in strain names.

Running xenoGI
--------------

xenoGI is a command line program that sometimes can take a while to run. If you are working on a remote machine, it may be useful to run xenoGI from within ``screen``, which is available on most linux distributions. For this tutorial, ``screen`` shouldn't be necessary because everything runs in a few minutes. But if you move on to larger datasets it might be helpful.

What screen does is provide a command line which you can "detach". You can then logout of the machine, and your process will keep running. When you log back in, you can retrieve it.

Parsing the gbff files
~~~~~~~~~~~~~~~~~~~~~~

The very first thing we'll do is have xenoGI run through these genbank files, and extract the protein annotations that we'll be using::

  python3 ../xenoGI/xenoGI-runner.py params.py parseGenbank

This should take 10-15 seconds.

Note that we wrote ``python3`` above, but on some systems you may want to write simply ``python``. Just be sure that this is calling the correct version of python, with the various necessary python packages. If you are using a pip-installed verion of xenoGI, then your command would look like this::

  xenoGI params.py parseGenbank

(You can make the equivalent adjustment for the commands to follow).

See README for a description of what fields are kept for each gene.

Running blast
~~~~~~~~~~~~~

Next we do an all vs. all protein blast::

  python3 ../xenoGI/xenoGI-runner.py params.py runBlast

This will take several minutes. For the steps below we will also try to give you a sense how long it should take on the tutorial data set. Note that speed may vary somewhat on your setup, but these numbers should give you a rough idea. If you subsequently do this on a larger data set of your own, of course it will take longer.

Calculating scores
~~~~~~~~~~~~~~~~~~

And then we calculate various types of scores::

  python3 ../xenoGI/xenoGI-runner.py params.py calcScores

The scores calulated in include raw similarity scores, as well as two types of synteny score. The core synteny score measures synteny in a large genomic neighborhood. The regular synteny score represent synteny in a smaller neighborhood. All three types of score can range from 0-1.
  
This step should take about 30 seconds.
  
Determining the species tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this step we'll determine the species tree for the strains we're looking at.

This step requires that the user specify an outgroup to root the species tree. In the enteric data set we're using, C_rodentium is the outgroup. Before we run the step, we need to specify the outgroup. In the 'Making species trees' section of ``params.py``, there is a parameter ``outGroup`` which has been commented out. Uncomment this (delete the hash) and set it so it reads::

  outGroup = 'C_rodentium'

Then run like so::

  python3 ../xenoGI/xenoGI-runner.py params.py makeSpeciesTree

This should take a minute or so, and will produce a newick file called ``enterics.tre``. Note that xenoGI doesn't make use of branch lengths, and the newick file produced here does not contain them.

For your reference, here's an ascii drawing of the resulting tree, with internal nodes labelled::

         ________________ E_coli_K12
    ____|
   |    |s1    __________ S_bongori
   |    |_____|
  _|          |s2   _____ S_enterica_LT2
   |s0        |____|s3
   |               |_____ S_enterica_AZ
   |
   |_____________________ C_rodentium


What you would do if you already knew the species tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When working on your own data, if you already know the tree, then you can enter it directly and skip ``makeSpeciesTree``. The species tree needs to be in newick format. It should have named internal nodes, and does not need to have branch lengths (if it has them, they will be ignored). The parameter ``speciesTreeFN`` in ``params.py`` gives the file name for this tree.

If using your own tree, make sure that the names in it match those being used by xenoGI (e.g. the names provided in ``ncbiHumanMap.txt``).

For reference, here's the newick string for the tree reconstructed above::

  ((E_coli_K12,(S_bongori,(S_enterica_LT2,S_enterica_AZ)s3)s2)s1,C_rodentium)s0;

The script ``prepareSpeciesTree.py`` in the ``misc/`` folder may be useful in preparing trees in the right format. (e.g. if you had a tree that was output from a phylogenetic reconstruction program). It will take an unrooted tree in newick and root it, as well as naming the internal nodes. It is described further in the documentation inside ``misc/``.
  
Creating gene families
~~~~~~~~~~~~~~~~~~~~~~

xenoGI does its most detailed reconstruction within a focal clade, leaving one or more species as outgroups. Such outgroups help us to better recognize core genes given the possibility of deletion in some lineages. One parameter we must set is the root of the focal clade. Once again, edit the ``params.py`` file. The line defining the ``rootFocalClade`` should be as follows::

  rootFocalClade = 's2'

This says that the focal clade will be defined by the internal node ``s2``, and corresponds to the Salmonella genus. ``C_rodentium`` and ``E_coli_K12`` will be outgroups.

We will now create gene families like so::

  python3 ../xenoGI/xenoGI-runner.py params.py makeFamilies

This will take several minutes on the tutorial data set.

``makeFamilies`` goes through several steps, ultimately producing a set of "origin" families which encompass genes with a common origin. In general, the possible types of origin are core gene or xeno hgt (horiztonal transfer from outside the clade). An origin family resulting from xeno hgt would include all the genes that are descended from one ancestor gene that arrived via hgt. A core origin family would include all genes descended from a core genea in the common ancestor strain. Origin families are broken up into locus families which consist of genes that occur in the same syntenic region. An origin family will consist of one or more locus families.

The process of making families involves, among other things, creating a gene trees and reconciling them against the species tree. For this we use the DTLOR reconciliation model. For more information see the README.

Creating locus islands
~~~~~~~~~~~~~~~~~~~~~~

Next create locus islands::
  
  python3 ../xenoGI/xenoGI-runner.py params.py makeIslands

This step involves grouping locus families that have a common origin into locus islands. It will likely take 1-2 minutes.

Refinement step
~~~~~~~~~~~~~~~

Finally we refine families and remake islands::

  python3 ../xenoGI/xenoGI-runner.py params.py refine

This will also take 1-2 minutes. In the refinement step, xenoGI goes back and looks at cases where there are multiple most-parsimonious reconciliations. In the previous ``makeFamilies`` step, one of these was chosen arbitrarily. Now xenoGI considers all of the possibilities, and determines which of these is optimal by examining nearby gene families. (On the logic that since these will often have a common origin, it makes sense to chose the most-parsimonious reconciliation the corresponds best to them.)

Creating output files
~~~~~~~~~~~~~~~~~~~~~

We can now create a set of output files which we'll use in subsequent analysis::

  python3 ../xenoGI/xenoGI-runner.py params.py printAnalysis

This step is very quick, taking just a few seconds on this data set.

Analysis
--------

Examining the genes files
~~~~~~~~~~~~~~~~~~~~~~~~~

The above command creates a subdirectory called analysis. Inside it you should find a set of files beginning with "genes", as well as ``islandsSummary.txt`` and ``islands.tsv``.

The genes files contain all the genes in a strain laid out in the order they occur on the contigs (the first line of each specifies what the columns are). Let's start out by looking at a known pathogenicity island, Salmonella Pathogenicity Island 1 (SPI1). This island is known to be present in all three Salmonella strains, S_enterica_LT2, S_enterica_AZ, and S_bongori. In S_enterica_LT2 it is known to extend from STM2865 to STM2900. Let's take a look using a text viewer. From within the analysis directory (``xgiTutorial/enterics/analysis/``) type::

  less -S genes-S_enterica_LT2.tsv

The ``-S`` tells the text viewer less not to wrap lines, which makes it a little easier to read. You may want to maximize your window, or make it wider so that more of each line displays. At the right of each line is included a description of each gene.

FYI, when you want to exit ``less``, type ``q``.

You can now search within ``less`` by typing forward slash (``/``) and entering the terms you want to search with. Here let's search using locus tag STM2865 which is at the beginning of SPI1.

Here's a truncated bit of what you should see::

  21087_S_enterica_LT2-STM2863  C       OSSS    3225    2724    3166    3225    s0      sitC - iron ABC transporter
  21088_S_enterica_LT2-STM2864  C       OSSS    3224    2723    3165    3224    s0      sitD - iron ABC transporter
  21089_S_enterica_LT2-STM2865  X       OS      3642    4170    4996    5069    s2      avrA - putative inner membr
  21090_S_enterica_LT2-STM2866  X       OSS     3642    3228    3788    3857    s2      sprB - transcriptional regu
  21091_S_enterica_LT2-STM2867  X       OSS     3642    3053    3578    3645    s2      hilC - AraC family transcri
  21092_S_enterica_LT2-STM2868  X       OSS     3642    3317    3901    3972    s2      type III secretion system e
  21093_S_enterica_LT2-STM2869  X       OSS     3642    3316    3900    3971    s2      orgA - invasion protein Org
  21094_S_enterica_LT2-STM2870  X       OSS     3642    3315    3899    3970    s2      putative inner membrane pro
  21095_S_enterica_LT2-STM2871  X       OSS     3642    3209    3759    3828    s2      prgK - EscJ/YscJ/HrcJ famil
  21096_S_enterica_LT2-STM2872  X       OSS     3642    3314    3898    3969    s2      prgJ - type III secretion s
  21097_S_enterica_LT2-STM2873  X       OSS     3642    3313    3897    3968    s2      prgI - EscF/YscF/HrpA famil
  21098_S_enterica_LT2-STM2874  X       OSS     3642    3312    3896    3967    s2      prgH - type III secretion s
  21099_S_enterica_LT2-STM2875  X       OSS     3642    3051    3576    3642    s2      hilD - AraC family transcri
  21100_S_enterica_LT2-STM2876  X       OSS     3642    3248    3816    3885    s2      hilA - transcriptional regu
  21101_S_enterica_LT2-STM2877  X       OSS     3642    3188    3737    3806    s2      iagB - invasion protein Iag
  21102_S_enterica_LT2-STM2878  X       OSS     3642    3311    3895    3966    s2      sptP - pathogenicity island
  21103_S_enterica_LT2-STM2879  X       OSS     3642    3310    3894    3965    s2      sicP - chaperone protein Si
  21104_S_enterica_LT2-STM2880  X       OS      4778    3941    4705    4778    s3      putative cytoplasmic protei
  21105_S_enterica_LT2-STM2881  X       OSS     3642    3160    3701    3770    s2      iacP - putative acyl carrie
  21106_S_enterica_LT2-STM2882  X       OSS     3642    3309    3893    3964    s2      sipA - pathogenicity island
  21107_S_enterica_LT2-STM2883  X       OSS     3642    3308    3892    3963    s2      sipD - cell invasion protei
  21108_S_enterica_LT2-STM2884  X       OSS     3642    3307    3891    3962    s2      sipC - pathogenicity island
  21109_S_enterica_LT2-STM2885  X       OSS     3642    3306    3890    3961    s2      sipB - pathogenicity island
  21110_S_enterica_LT2-STM2886  X       OSS     3642    3187    3736    3805    s2      sicA - CesD/SycD/LcrH famil
  21111_S_enterica_LT2-STM2887  X       OSS     3642    3126    3657    3726    s2      spaS - EscU/YscU/HrcU famil
  21112_S_enterica_LT2-STM2888  X       OSS     3642    3208    3758    3827    s2      spaR - EscT/YscT/HrcT famil
  21113_S_enterica_LT2-STM2889  X       OSS     3642    3207    3757    3826    s2      spaQ - EscS/YscS/HrcS famil
  21114_S_enterica_LT2-STM2890  X       OSS     3642    3125    3656    3725    s2      spaP - EscR/YscR/HrcR famil
  21115_S_enterica_LT2-STM2891  X       OSS     3642    3305    3889    3960    s2      spaO - type III secretion s
  21116_S_enterica_LT2-STM2892  X       OSS     3642    3304    3888    3959    s2      invJ - antigen presentation
  21117_S_enterica_LT2-STM2893  X       OSS     3642    3303    3887    3958    s2      invI - type III secretion s
  21118_S_enterica_LT2-STM2894  X       OSS     3642    3058    3583    3651    s2      invC - EscN/YscN/HrcN famil
  21119_S_enterica_LT2-STM2895  X       OSS     3642    3302    3886    3957    s2      invB - type III secretion s
  21120_S_enterica_LT2-STM2896  X       OSS     3642    3124    3655    3724    s2      invA - EscV/YscV/HrcV famil
  21121_S_enterica_LT2-STM2897  X       OSS     3642    3301    3885    3956    s2      invE - SepL/TyeA/HrpJ famil
  21122_S_enterica_LT2-STM2898  X       OSS     3642    3206    3756    3825    s2      invG - EscC/YscC/HrcC famil
  21123_S_enterica_LT2-STM2899  X       OSS     3642    3300    3884    3955    s2      invF - invasion protein
  21124_S_enterica_LT2-STM2900  X       OSS     3642    3299    3883    3954    s2      invH - invasion lipoprotein
  21125_S_enterica_LT2-STM2901  X       O       4588    3864    4603    4676    S_enterica_LT2  hypothetical protei
  21126_S_enterica_LT2-STM2902  X       O       4588    3802    4515    4588    S_enterica_LT2  putative cytoplasmi

The first column consists of genes listed by their xenoGI name (the locus tag is the last part of this). xenoGI has identified a locus island that corresponds to SPI1. The number for this locus island is given in column 4, and is 3646 here. (It is possible that the numbering will be different on your machine). This locus island extends from 21089_S_enterica_LT2-STM2865 to 21124_S_enterica_LT2-STM2900 as expected. Note that in the display above, we've included a few genes on either end of the locus island.

As discussed in the README, a locus island represents a set of gene families with a common origin. In this case, it corresponds to a genomic island which is inferred to have inserted on the branch leading to s2 (the branch inserted on is given in the 8th column).

Every gene in a particular clade is either a core gene, or arose by xeno horizontal transfer (horizontal transfer from outside the clade). One of the goals of xenoGI is to determine this origin for each gene. The second column in the genes file contains this information. C stands for core, and X for xeno horizontal transfer. Note that for SPI1, all the genes are marked X.

The third column contains a gene history string. Taking the gene 21124_S_enterica_LT2-STM290 for example (invH) the string is OSS. This reflects the history of the gene after insertion, as reconstructed by the DTLOR reconciliation. O stands for origin (in this case the xeno hgt event). And S stands for co-speciation--what happens when a speciation event occurs and both descendent lineages inherit a gene. invH is inferred to have inserted on branch s2. It then underwent co-speciation events at node s2 and node s3. Other possible characters that could appear in the gene history string are  D, duplication; T, transfer (within the species tree); R, rearrangement (from one syntenic region to another).

As we noted, the 4th column gives the locus island. The 5th gives the initial family number, the 6th the origin family number, and the 7th the locus family number. We'll use some of these in the examples below.

Quit out of ``less`` by typing ``q``.

A second pathogenicity island in Salmonella, SPI2 is known to have two parts with different evolutionary origins. The type III secretion system (t3ss) is shared by Salmonella enterica strains, but is lacking outside that group. On our enterics tree, this means it inserted on the s3 branch. There is also a portion of SPI2 that is called the tetrathionate reductase gene cluster (trgc). This portion is present in other species in the Salmonella genus. On our enterics tree it inserted on the s2 branch. The following locus tags define the beginning and end of these regions in SPI2 in S_enterica_AZ.

==== ========== ==========
 \     From         To
==== ========== ==========
t3ss SARI_01560 SARI_01590
trgc SARI_01591 SARI_01600
==== ========== ==========

You can search for these as we did above, and see what xenoGI says about the origins of these genes::

  less -S genes-S_enterica_AZ.tsv

Examining island summary filess
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's now take a look at a second file::

  less islandsSummary.txt

This file provides a human readable listing of locus islands, organized by the branch where they inserted. If you search for "LocusIsland 3646" it will bring you to the entry for the SPI1 island. Each entry has two parts. First is a listing of families, written out by row. Then below that is a listing of the genes that includes the description of the gene.

This file is especially useful if you are browsing for interesting novel islands.

Note that there is a tab delimited version of this information contained in the file ``islands.tsv`` (which will be more useful if you want to read it in to some subsequent analysis program).

Interactive analysis
~~~~~~~~~~~~~~~~~~~~

It is possible to get additional information using the interactive analysis mode.

Let's say you have a gene and are wanting to learn about the evolution of the family it belongs to. For example, maybe you are interested in SARI_01595 which is part of the tetrathionate reductase gene cluster. To proceed, you need to know the origin family number for this gene.

One way to find out is to look in the genes files, as described above. Another way involves using the ``findGene`` function in interactive analysis.

From a terminal prompt in the main enterics directory (``xgiTutorial/enterics``) type::

  python3 ../xenoGI/xenoGI-runner.py params.py interactiveAnalysis

Once the python promp comes up, type::

  findGene("SARI_01595")

(Note that tab completion works in interactive mode, at least on Linux).

The output should look like this::

  <gene:1534_S_enterica_AZ-SARI_01595 locIsl:2006 ifam:3506 ofam:4143 locFam:4216 hypothetical protein>

``findGene`` searches all the information associated with a gene. So you can potentially give it not only a locus tag, but also common name, protein ID etc.
  
``printFam``
^^^^^^^^^^^^
We can see from this that the origin family is 4143. (It is possible that on your machine the numbers will be different. If so, substitute the number you got for 4143 below).

We can now type the following at the python prompt::

  printFam(originFamiliesO,4143)
  
This produces the following output::

    Family 4143
        LocusFamily 4216 s2 4189 root_b 1534_S_enterica_AZ-SARI_01595 19655_S_enterica_LT2-STM1387 14955_S_bongori-A464_1417

        Source family 3506


    Matrix of raw similarity scores [0,1] between genes in the family
                                    | 1534_S_enterica_AZ-SARI_01595 | 19655_S_enterica_LT2-STM1387 | 14955_S_bongori-A464_1417
      1534_S_enterica_AZ-SARI_01595 | 1.000                         | 0.944                        | 0.896
      19655_S_enterica_LT2-STM1387  | 0.944                         | 1.000                        | 0.919
      14955_S_bongori-A464_1417     | 0.896                         | 0.919                        | 1.000


    Matrix of core synteny scores [0,1] between genes in the family
                                    | 1534_S_enterica_AZ-SARI_01595 | 19655_S_enterica_LT2-STM1387 | 14955_S_bongori-A464_1417
      1534_S_enterica_AZ-SARI_01595 | 1.000                         | 1.000                        | 1.000
      19655_S_enterica_LT2-STM1387  | 1.000                         | 1.000                        | 1.000
      14955_S_bongori-A464_1417     | 1.000                         | 1.000                        | 1.000


    Matrix of synteny scores [0,1] between genes in the family
                                    | 1534_S_enterica_AZ-SARI_01595 | 19655_S_enterica_LT2-STM1387 | 14955_S_bongori-A464_1417
      1534_S_enterica_AZ-SARI_01595 | 1.000                         | 0.972                        | 0.940
      19655_S_enterica_LT2-STM1387  | 0.972                         | 1.000                        | 0.958
      14955_S_bongori-A464_1417     | 0.940                         | 0.958                        | 1.000


    Printing all scores with non-family members
      Inside fam                    | Outside fam                   | Raw   | Syn   | CoreSyn
      ----------                    | -----------                   | ---   | ---   | -------
      19655_S_enterica_LT2-STM1387  | 12484_C_rodentium-ROD_RS20160 | 0.398 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 12484_C_rodentium-ROD_RS20160 | 0.396 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 7970_E_coli_K12-b3669         | 0.396 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 7970_E_coli_K12-b3669         | 0.395 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 19659_S_enterica_LT2-STM1391  | 0.394 | 1.000 | 0.950
      14955_S_bongori-A464_1417     | 7970_E_coli_K12-b3669         | 0.394 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 12484_C_rodentium-ROD_RS20160 | 0.394 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 5657_E_coli_K12-b1221         | 0.391 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 21014_S_enterica_LT2-STM2785  | 0.390 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 21992_S_enterica_LT2-STM3790  | 0.389 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 17399_S_bongori-A464_3863     | 0.389 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 1529_S_enterica_AZ-SARI_01590 | 0.389 | 0.638 | 0.950
      14955_S_bongori-A464_1417     | 182_S_enterica_AZ-SARI_00190  | 0.389 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 19659_S_enterica_LT2-STM1391  | 0.389 | 0.307 | 0.950
      14955_S_bongori-A464_1417     | 10345_C_rodentium-ROD_RS08835 | 0.389 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 20203_S_enterica_LT2-STM1947  | 0.389 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 21992_S_enterica_LT2-STM3790  | 0.388 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 20026_S_enterica_LT2-STM1767  | 0.388 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 19659_S_enterica_LT2-STM1391  | 0.387 | 0.928 | 0.950
      14955_S_bongori-A464_1417     | 10520_C_rodentium-ROD_RS09755 | 0.387 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 15401_S_bongori-A464_1864     | 0.387 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 946_S_enterica_AZ-SARI_00990  | 0.387 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 21992_S_enterica_LT2-STM3790  | 0.387 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 1136_S_enterica_AZ-SARI_01186 | 0.387 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 17399_S_bongori-A464_3863     | 0.385 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 15634_S_bongori-A464_2097     | 0.385 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 17399_S_bongori-A464_3863     | 0.385 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 1529_S_enterica_AZ-SARI_01590 | 0.384 | 1.000 | 0.950
      14955_S_bongori-A464_1417     | 1529_S_enterica_AZ-SARI_01590 | 0.384 | 0.000 | 0.950
      1534_S_enterica_AZ-SARI_01595 | 5657_E_coli_K12-b1221         | 0.381 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 15401_S_bongori-A464_1864     | 0.381 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 20026_S_enterica_LT2-STM1767  | 0.381 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 1136_S_enterica_AZ-SARI_01186 | 0.381 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 5657_E_coli_K12-b1221         | 0.381 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 15401_S_bongori-A464_1864     | 0.381 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 10345_C_rodentium-ROD_RS08835 | 0.381 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 6328_E_coli_K12-b1914         | 0.381 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 20026_S_enterica_LT2-STM1767  | 0.380 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 1136_S_enterica_AZ-SARI_01186 | 0.380 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 20203_S_enterica_LT2-STM1947  | 0.379 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 10520_C_rodentium-ROD_RS09755 | 0.379 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 10345_C_rodentium-ROD_RS08835 | 0.378 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 10520_C_rodentium-ROD_RS09755 | 0.377 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 946_S_enterica_AZ-SARI_00990  | 0.377 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 20203_S_enterica_LT2-STM1947  | 0.376 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 7352_E_coli_K12-b3025         | 0.376 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 16744_S_bongori-A464_3208     | 0.376 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 15634_S_bongori-A464_2097     | 0.375 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 946_S_enterica_AZ-SARI_00990  | 0.375 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 7352_E_coli_K12-b3025         | 0.374 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 6328_E_coli_K12-b1914         | 0.373 | 0.000 | 0.000
      19655_S_enterica_LT2-STM1387  | 3698_S_enterica_AZ-SARI_03859 | 0.364 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 3698_S_enterica_AZ-SARI_03859 | 0.363 | 0.000 | 0.000
      1534_S_enterica_AZ-SARI_01595 | 3698_S_enterica_AZ-SARI_03859 | 0.362 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 10860_C_rodentium-ROD_RS11510 | 0.360 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 18847_S_enterica_LT2-STM0549  | 0.354 | 0.000 | 0.000
      14955_S_bongori-A464_1417     | 6749_E_coli_K12-b2369         | 0.354 | 0.000 | 0.000


    Gene tree
    ((1534,19655)g0,14955)root

    Gene tree annotated with reconciliation [branch events | node events]
    ((1534[|S_enterica_AZ],19655[|S_enterica_LT2])g0[|S],14955[|S_bongori])root[O|S]

    Reconciliation of gene tree onto species tree
    - root
      O (root b) --> (s2 b) synReg:4189
      S (root n) --> (s2 n)
       - g0
         S (g0 n) --> (s3 n)
          - 1534 [S_enterica_AZ]
          - 19655 [S_enterica_LT2]
       - 14955 [S_bongori]

Note that if you want to save this output directly to a file you can do like this::

    printFam(originFamiliesO,4143,open("ofam4143.txt","w"))

The third argument is optional, and is an open file handle. Doing this can be useful if you have a large family, and you want to view it without lines wrapping. (e.g. with ``less -S``).

Let's now go though the various parts of this output.

The family we've just printed is an origin family. Origin families represent the more refined stage of family analysis, and are what users are most likely to be interested in. An origin family has a gene tree associated with it, and also a reconciliation that places that gene tree onto the species tree. At the base of this reconciliation is an origin event. In this case, it is a xeno hgt event. (The other alternative is if a family is a core gene family).

The first line of the output gives the family number.

Next come some lines printing out the locus families that are part of this origin family. A locus family represents the genes in a family which occur in a single syntenic region. Every family has at least one locus family, but may have more. 

In this case there is only one locus family, number 4216. This locus family originated on branch s2 of the species tree. It it found in syntenic region 4189, and it's origin point in the gene tree is the root branch of that tree. The remainder of this line consists of a listing of the genes in this locus family. In this case there are 3, one in each of the 2 S. enterica strains, and 1 in the S. bongori strain.

The next element of the output is the source family. Because we're looking at an origin family, the source family for that will be what we call an initial family--initial family 3506 in this case. (This bit of information would be useful if you wanted to go back to the initial family and look at how it was split up into origin families).

The next elements are 3 score matricies, showing the raw, core synteny, and regular synteny scores for the genes in this origin family. All of these scores take on values ranging from 0 to 1.

The raw score is a sequence similarity score. In this case, all 3 genes are fairly similar.

The core synteny score reflects synteny as defined relative to core genes (large scale or long distance synteny). In this case, we can see that all 3 genes are in the same location given the very high scores.

The regular synteny score represents more fine grained synteny, looking at a neighborhood of 20 genes around each family member. These synteny scores are also high in this case.

The high synteny between all family members is the reason that xenoGI made only a single locus family in this origin family.

The next element of the output is a printout of scores between family members, and non-family members. (The non-family members represent all genes that have significant blast hits vs. family members.) You might be interested in this if you suspected that there were some genes left out of the family that should have been included.

In this case, all non family members are very different in terms of their sequences (low raw scores). Most of them also reside in different syntenic regions. There is nothing in this list that looks like a gene which should have been included in this family.

The next elements of the output are a gene tree in newick format, and also an "annotated" version of the gene tree. One way to view these is to cut and paste the newick string into a file, add a semicolon at the end, and then view this with a tree viewer such as FigTree.

If you use FigTree in this way, it imports the annotation (called "label" by default) and then lets you display it on the nodes. For example, at the root of the species tree, there is this annotation::

  root[O|S]

Inside the bracket, we have two elements separated by a ``|``. The left one represents events that happened on the branch leading to root, and the right one represents events that happened on the node. Here, we have an origin event on the the branch leading to the root (in this case, a xeno hgt event). At the root node of the gene tree, we have a co-speciation event, where the species tree diverges, and each descendent lineage inherits a copy of the gene. On tips of the gene tree, the node part of this (part on the right) will simply give the strain name of the strain where the tip gene is found.

The final element of the output is a text representation of the reconciliation. This representation is organized according to the gene tree. So it basically goes through the gene tree, and specifies events occuring on gene tree branches and nodes, and the placements onto the species tree.

We begin with the root of the gene tree. There is a listing of events. An O (origin) event occurred on the root branch of the gene tree and the s2 branch on the species tree. Because s2 is in the focal clade, an O event here really represents xeno hgt. The text representing this event also tells us that the insertion occurred at syntenic region 4189. The second event listed is a S event (cospeciation). This involves the placement of the root node of the gene tree on the s2 node of the species tree. So there was a cospeciation at s2 where the gene was interited in each of the two species tree lineages descending from s2.

Listed below this is what happened to the two children of the gene tree root node, g0 and gene 14955. Gene 14955 is a tip on the gene tree, and is found in S_bongori. At g0, there was another S event (cospeciation). g0 is placed on the species node s3. There is a cospeciation event there where the two descendent branches of the gene tree, genes 1534 and 19655, are inherited in S_enterica_AZ and S_enterica_LT2 respectively.

``printLocusIsland``
^^^^^^^^^^^^^^^^^^^^

Sometimes you might be interested in looking at a particular locus island, and seeing it in each of the strains where it occurs. One way to do this is to look through all the genes files for those strains (as described above).

However, interactive analysis provides a convenient way of printing a locus island in all the strains where it occurs.

For example, the gene SARI_01595 that we were interested in above is part of locus island 2006. Let's view that locus island.

At the python prompt (which you got by running the interactiveAnalysis command) type the following::

  printLocusIsland(2006,20)

This will print locus island 2006, showing 20 genes surrounding (10 in either direction). Those genes that are part of locus island 2006 are indicated with a star. The columns are the same as what is in the genes files, as described above. Also included are the genomic coordinates of the island and the region::

  LocusIsland: 2006
  mrca: s2
  In S_bongori
    Coordinates of locus island CP006608.1:1401420-1409685
    Coordinates of region shown CP006608.1:1394919-1419664
    geneName                    | orig | geneHist | locIsl | ifam | ofam | locFam | lfMrca    | descrip
      14943_S_bongori-A464_1405 | C    | OSSS     | 2008   | 1573 | 1957 | 2008   | s0        | Iron-sulfur cluster assembly protein SufD
      14944_S_bongori-A464_1406 | C    | OSSS     | 905    | 528  | 864  | 905    | s0        | Cysteine desulfurase subunit
      14945_S_bongori-A464_1407 | C    | OSSS     | 1228   | 834  | 1186 | 1228   | s0        | Sulfur acceptor protein SufE for iron-sulfurcluster assembly
      14946_S_bongori-A464_1408 | C    | OSSS     | 772    | 402  | 731  | 772    | s0        | LD-transpeptidase YnhG
      14947_S_bongori-A464_1409 | C    | OSDS     | 414    | 88   | 381  | 414    | s0        | major outer membrane lipoprotein
      14948_S_bongori-A464_1410 | X    |          | 9126   | 8204 | 9053 | 9126   | S_bongori | hypothetical protein
      14949_S_bongori-A464_1411 | C    | OSSS     | 1184   | 791  | 1142 | 1184   | s0        | Pyruvate kinase
    * 14950_S_bongori-A464_1412 | X    | OS       | 2006   | 3990 | 4761 | 4834   | s2        | Putative amino acid permease
    * 14951_S_bongori-A464_1413 | X    | OS       | 2006   | 3374 | 3960 | 4032   | s2        | Tetrathionate reductase subunit A
    * 14952_S_bongori-A464_1414 | X    | OS       | 2006   | 3373 | 3959 | 4031   | s2        | Tetrathionate reductase subunit C
    * 14953_S_bongori-A464_1415 | X    | OS       | 2006   | 3037 | 3554 | 3620   | s2        | Tetrathionate reductase subunit B
    * 14954_S_bongori-A464_1416 | X    | OS       | 2006   | 3372 | 3958 | 4030   | s2        | Tetrathionate reductase sensory transductionhistidine kinase
    * 14955_S_bongori-A464_1417 | X    | OS       | 2006   | 3506 | 4143 | 4216   | s2        | Tetrathionate reductase two-component responseregulator
    * 14956_S_bongori-A464_1418 | X    | OS       | 2006   | 1572 | 1955 | 2006   | s2        | hypothetical protein
      14957_S_bongori-A464_1419 | X    |          | 3590   | 4326 | 5175 | 5248   | S_bongori | Transcriptional regulatory protein
      14958_S_bongori-A464_1420 | X    | O        | 3590   | 3011 | 3525 | 3590   | S_bongori | Alcohol dehydrogenase
      14959_S_bongori-A464_1421 | X    | ODS      | 3397   | 2860 | 3335 | 3397   | s2        | hypothetical protein
      14960_S_bongori-A464_1422 | X    | OS       | 3397   | 3116 | 3646 | 3715   | s2        | Transcriptional regulator MerR familyassociated with photolyase
      14961_S_bongori-A464_1423 | X    |          | 1557   | 5993 | 6842 | 6915   | S_bongori | hypothetical protein
      14962_S_bongori-A464_1424 | X    |          | 1557   | 5995 | 6844 | 6917   | S_bongori | Uncharacterized protein ImpA
      14963_S_bongori-A464_1425 | X    |          | 1557   | 8205 | 9054 | 9127   | S_bongori | IcmF-related protein
  In S_enterica_LT2
    Coordinates of locus island NC_003197.2:1466345-1474023
    Coordinates of region shown NC_003197.2:1459047-1483078
    geneName                       | orig | geneHist | locIsl | ifam | ofam | locFam | lfMrca         | descrip
      19644_S_enterica_LT2-STM1376 | C    | OSDSS    | 414    | 88   | 381  | 414    | s0             | lppB - hypothetical protein
      19645_S_enterica_LT2-STM1377 | C    | OSDS     | 414    | 88   | 381  | 414    | s0             | lpp - murein lipoprotein
      19646_S_enterica_LT2-STM1378 | C    | OSSSS    | 1184   | 791  | 1142 | 1184   | s0             | pykF - pyruvate kinase
      19647_S_enterica_LT2-STM1379 | X    |          | 5315   | 4898 | 5747 | 5820   | S_enterica_LT2 | orf48 - putative amino acid permease
      19648_S_enterica_LT2-STM1380 | X    |          | 5315   | 8829 | 9678 | 9751   | S_enterica_LT2 | orf32 - hydrolase
      19649_S_enterica_LT2-STM1381 | X    |          | 5315   | 5620 | 6469 | 6542   | S_enterica_LT2 | orf245 - hypothetical protein
      19650_S_enterica_LT2-STM1382 | X    |          | 5315   | 4393 | 5242 | 5315   | S_enterica_LT2 | orf408 - hypothetical protein
    * 19651_S_enterica_LT2-STM1383 | X    | OSS      | 2006   | 3374 | 3960 | 4032   | s2             | ttrA - tetrathionate reductase subunit A
    * 19652_S_enterica_LT2-STM1384 | X    | OSS      | 2006   | 3373 | 3959 | 4031   | s2             | ttrC - tetrathionate reductase subunit C
    * 19653_S_enterica_LT2-STM1385 | X    | OSS      | 2006   | 3037 | 3554 | 3620   | s2             | ttrB - tetrathionate reductase complex, subunit B
    * 19654_S_enterica_LT2-STM1386 | X    | OSS      | 2006   | 3372 | 3958 | 4030   | s2             | ttrS - tetrathionate reductase complex: sensory transduction histidine kinase
    * 19655_S_enterica_LT2-STM1387 | X    | OSS      | 2006   | 3506 | 4143 | 4216   | s2             | ttrR - DNA-binding response regulator
    * 19656_S_enterica_LT2-STM1388 | X    | OSS      | 2006   | 1572 | 1955 | 2006   | s2             | orf70 - hypothetical protein
      19657_S_enterica_LT2-STM1389 | X    | OD       | 3397   | 2860 | 3335 | 3397   | s2             | orf319 - hypothetical protein
      19658_S_enterica_LT2-STM1390 | X    | OSS      | 3397   | 3116 | 3646 | 3715   | s2             | orf242 - helix-turn-helix-type transcriptional regulator
      19659_S_enterica_LT2-STM1391 | X    | OS       | 4349   | 4210 | 5058 | 5131   | s3             | ssrB - DNA-binding response regulator
      19660_S_enterica_LT2-STM1392 | X    | OS       | 4349   | 3989 | 4760 | 4833   | s3             | ssrA - hybrid sensor histidine kinase/response regulator
      19661_S_enterica_LT2-STM1393 | X    | OS       | 4349   | 3988 | 4759 | 4832   | s3             | ssaB - pathogenicity island chaperone protein SpiC
      19662_S_enterica_LT2-STM1394 | X    | OS       | 4349   | 3768 | 4468 | 4541   | s3             | ssaC - EscC/YscC/HrcC family type III secretion system outer membrane ring protein
      19663_S_enterica_LT2-STM1395 | X    | OS       | 4349   | 3878 | 4619 | 4692   | s3             | ssaD - EscD/YscD/HrpQ family type III secretion system inner membrane ring protein
  In S_enterica_AZ
    Coordinates of locus island CP000880.1:1542713-1550949
    Coordinates of region shown CP000880.1:1534729-1558172
    geneName                        | orig | geneHist | locIsl | ifam | ofam | locFam | lfMrca        | descrip
      1526_S_enterica_AZ-SARI_01587 | X    | OS       | 4349   | 3768 | 4468 | 4541   | s3            | hypothetical protein
      1527_S_enterica_AZ-SARI_01588 | X    | OS       | 4349   | 3988 | 4759 | 4832   | s3            | hypothetical protein
      1528_S_enterica_AZ-SARI_01589 | X    | OS       | 4349   | 3989 | 4760 | 4833   | s3            | hypothetical protein
      1529_S_enterica_AZ-SARI_01590 | X    | OS       | 4349   | 4210 | 5058 | 5131   | s3            | hypothetical protein
      1530_S_enterica_AZ-SARI_01591 | X    | OSS      | 3397   | 3116 | 3646 | 3715   | s2            | hypothetical protein
      1531_S_enterica_AZ-SARI_01592 | X    | ODS      | 3397   | 2860 | 3335 | 3397   | s2            | hypothetical protein
      1532_S_enterica_AZ-SARI_01593 | X    |          | 7309   | 6387 | 7236 | 7309   | S_enterica_AZ | hypothetical protein
    * 1533_S_enterica_AZ-SARI_01594 | X    | OSS      | 2006   | 1572 | 1955 | 2006   | s2            | hypothetical protein
    * 1534_S_enterica_AZ-SARI_01595 | X    | OSS      | 2006   | 3506 | 4143 | 4216   | s2            | hypothetical protein
    * 1535_S_enterica_AZ-SARI_01596 | X    | OSS      | 2006   | 3372 | 3958 | 4030   | s2            | hypothetical protein
    * 1536_S_enterica_AZ-SARI_01597 | X    | OSS      | 2006   | 3037 | 3554 | 3620   | s2            | hypothetical protein
    * 1537_S_enterica_AZ-SARI_01598 | X    | OSS      | 2006   | 3373 | 3959 | 4031   | s2            | hypothetical protein
    * 1538_S_enterica_AZ-SARI_01599 | X    | OSS      | 2006   | 3374 | 3960 | 4032   | s2            | hypothetical protein
      1539_S_enterica_AZ-SARI_01601 | X    |          | 7310   | 6388 | 7237 | 7310   | S_enterica_AZ | hypothetical protein
    * 1540_S_enterica_AZ-SARI_01600 | X    | OS       | 2006   | 3990 | 4761 | 4834   | s2            | hypothetical protein
      1541_S_enterica_AZ-SARI_01602 | C    | OSSSS    | 1184   | 791  | 1142 | 1184   | s0            | hypothetical protein
      1542_S_enterica_AZ-SARI_01603 | C    | OSDSS    | 414    | 88   | 381  | 414    | s0            | hypothetical protein
      1543_S_enterica_AZ-SARI_01604 | C    | OSSSS    | 772    | 402  | 731  | 772    | s0            | hypothetical protein
      1544_S_enterica_AZ-SARI_01605 | C    | OSSSS    | 1228   | 834  | 1186 | 1228   | s0            | hypothetical protein
      1545_S_enterica_AZ-SARI_01606 | C    | OSSSS    | 905    | 528  | 864  | 905    | s0            | hypothetical protein
      1546_S_enterica_AZ-SARI_01607 | C    | OSSSS    | 2008   | 1573 | 1957 | 2008   | s0            | hypothetical protein
      1547_S_enterica_AZ-SARI_01608 | C    | OSSSS    | 3064   | 2580 | 3007 | 3064   | s0            | hypothetical protein

Locus island 2006 corresponds to the tetrathionate reductase gene cluster. In fact, several additional locus families (3397,3715) probably should have been included in locus island 2006. They likely all had a common origin. They reason xenoGI did not include them is that several strain specific genes have been inserted between them and the rest of the island (genes 14957_S_bongori-A464_1419, and 14958_S_bongori-A464_1420 in S_bongori, and 1532_S_enterica_AZ-SARI_01593 in S_enterica_AZ). This illustrates a limitation: xenoGI tries to group everythig with a common origin, but sometimes the evolutionary history makes it hard to do that.

Note that in the S_enterica species you can also see the nearby type III secretion system, which is island 4349 (not shown in its entirety).


Viewing in a genome browser
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can create bed files of the output and then view them in a genome browser. This makes possible a representation where locus islands are colorized. Create the beds::

  python3 ../xenoGI/xenoGI-runner.py params.py createIslandBed

This creates a ``bed`` subdirectory with a bed file for each strain.

Such files can be viewed with a variety of browsers. Here we'll give an example using the Ingegrated Genome Browser (IGB).

IGB can be downloaded here: https://bioviz.org/ .

In this example we're using version 9.1.6 (other versions will likely work as well).

We're going to view the S_enterica_LT2 genome. Therefore, move the ``S_enterica_LT2-island.bed`` file so that it's on the same machine where you're going to run IGB (ie if its on a remote machine, move it back to your local machine). Next go to the NCBI assembly page for this genome https://www.ncbi.nlm.nih.gov/assembly/GCF_000006945.2/ . Click on the link to the FTP directory for the refseq assembly (that happens to be the one we're using for LT2 in our example.) From the FTP directory, download GCF_000006945.2_ASM694v2_cds_from_genomic.fna.gz and GCF_000006945.2_ASM694v2_genomic.gff.gz.

Move these to the same location as the bed file. Then unzip them::

  gunzip GCF_000006945.2_ASM694v2_cds_from_genomic.fna.gz GCF_000006945.2_ASM694v2_genomic.gff.gz

Now start the IGB broswer. From the File menu, select "Open Genome from File". Select the fna file you just downloaded and unpacked. Once that has completed, again go to the File menu. This time select "Open File", and select ``S_enterica_LT2-island.bed``.

To get the genes to display, you may have to click "Load Data" in the upper right.

To make the islands show up in color, right click (or control click with a single button mouse) the blue control box for the bed track (it's on the left). Select the "Color By" option. Then when a window pops up asking you what to yous for coloring, select RGB and hit OK.

Now, let's go have a look at the SPI1 island. We can determine the coordinates for SPI1 from the printLocusIsland command in interactive mode. If you want, you can go back and do that now. Alternatively, here are the coordinates::

  NC_003197.2:3009904-3044839

You should paste this into the IGB coordinate window (upper left).

This will take you to the exact region of the locus island xenoGI found which corresponds to SPI1. You may want to zoom out a bit so you can see it in context.


