===============
xenoGI tutorial
===============


Introduction
------------

The purpose of this tutorial is to provide an example of how run xenoGI, and how to use its output. The tutorial is intended to be run on Linux or Mac. (But could probably work on Windows with some modifications).

We will assume some familiartity with the unix command line, but will otherwise try to work at a basic level.

For consistency in terms of file paths we will set up a tutorial directory and place a copy of the xenoGI repository, as well as the data set we're working on within it. To do this, start a terminal on the machine you're going to work on. In whatever location suits you, create the following directory::

  mkdir xgiTutorial

Then move into it::

  cd xgiTutorial

Next we're going to put a copy of the repository inside this directory. If you already have a copy, feel free to move it here (e.g. with the Mac gui, with the ``mv`` command etc). Alternatively you can clone a new copy. If you're on a Mac you could do this with a browser via the github website (https://github.com/ecbush/xenoGI). You can also do it at the command line (assuming you have git installed)::

  git clone https://github.com/ecbush/xenoGI


Getting the genomes
-------------------

We'll also create a directory for the data set we're going to work with::

  mkdir enterics
  cd enterics

xenoGI makes heavy use of synteny, and therefore requires genomes that are assembled to the scaffold level or better. We're going to work with a set of complete assemblies from enteric bacteria (this is a different set than the one distributed with the repository in the example/ directory)

Withing the enterics directory we'll create a subdirectory to hold the gbff files we're going to get::

  mkdir ncbi

