Code for detecting genomic island insertions in clades of bacteria.

## Requirements

- NCBI blast+
  We need a blastp executable in the path.

- Python 3, with

  . biopython
  http://biopython.org/
  For parsing genbank files.

  . parasail
  https://github.com/jeffdaily/parasail
  This is an optimized alignment library, used in calculating scores between proteins.


## How to use

An example/ directory is included in this repository.

We work with a set of species with known phylogenetic relationships. In the example, these species are: E. coli K12, E. albertii, E. fergusonii and S. bongori.

- The working directory must contain:

  . A parameter file. In the provided example/ directory this is called params.py.

  . A subdirectory of sequence files. In the example, this is called ncbi. Contained in this subdirectory will be genbank (gbff) files for the species.
     The parameter genbankFilePath in params.py has the path to these files.

  . A newick format tree representing the relationships of the strains. In the example this is called example.tre. Note that branch lengths are not used in xenoGI, and example.tre does not contain branch lengths. Also note that internal nodes should be given names in this tree. In the example.tre we label them i0, i1 etc. The parameter treeFN in params.py has the path to this tree file.

- Naming of genbank files.
  . The system needs a way to connect the sequence files to the names used in the tree.

   In the example, the sequence files have names corresponding to their assembly accession number from ncbi. We connect these to the human readable names in example.tre using a mapping given in the file ncbiHumanMap.txt. This file has two columns, the first giving the name of the genbank file, and the second giving the name for the species used in the tree file. Note that the species name should not contain any dashes ("-"). In params.py the parameter fileNameMapFN is set to point to this file.

   Another approach is to change the names of the sequence files to match what's in the tree. If you do this, then you should set fileNameMapFN = None in params.py. (This is not necessary in the example, which is already set to run the other way).

- Running the code
  
You run the code from within the working directory. To run the example, you would cd into the example/ directory, then:

python3 ../parseGenbank.py params.py
python3 ../runBlast.py params.py
python3 ../calcScores.py params.py
python3 ../xenoGI.py params.py

From another directory, this would simply be

python3 path-to-xenoGI-directory/parseGenbank.py params.py
python3 path-to-xenoGI-directory/runBlast.py params.py
python3 path-to-xenoGI-directory/calcScores.py params.py
python3 path-to-xenoGI-directory/xenoGI.py params.py


  . parseGenbank.py runs through the genbank files and produces input files that are used by subsequent code.
  . runBlast.py does an all vs. all protein blast of the genes in these strains. The number of processes it will run in parallel is specified by the numThreads parameter in the parameter file.
  . calcScores.py calculates similarity and synteny scores between genes in the strains. It is also (mostly) parallelized.
  . xenoGI.py identifies horizontal transfer events. It is partly parallelized.


- Subsequent analysis can also be run from the working directory

python3 -i path-to-xenoGI-directory/runAnalysis.py params.py

From within python, you can then run functions such as

  . printIslandsAtNode

     printIslandsAtNode('i0') # Prints all islands at node i0
     printIslandsAtNode('E_coli_K12') # Prints all islands specific to the E. coli K12 branch

  . findIsland 


    findIsland('gadA') # Find an island associated with a gene name or description

  . printIsland

    printIsland(3500,10) # Print out an island with its context in different species (first arg is island id, second is the number of genes to print to each side)

    printIsland prints the island in each strain where it's present. Its output includes the island and family numbers for each gene, an error score for the family of each gene, the most recent common ancestor (mrca) of the family, and a description of the gene. The error score is intended to indicate confidence in the correctness of the family. 0 means more confident, higher numbers less confident.

  . printFam

    printFam(3500) # Print scores within a particular gene family, and also with similar genes not in the family


- Visualization in a browser

  . We also include code to output the islands for each strain into a bed or gff file.

     python3 path-to-xenoGI-directory/misc/createIslandBed.py params.py 100

     python3 path-to-xenoGI-directory/misc/createIslandGffs.py params.py

  In the example, these will be created in a directory called bed/ or gff/ respectively.

  These can be visualized in a browser.

  . With the IGB browser (http://bioviz.org/igb/), the bed file option allows you to most easily display different islands in different colors.

  . To use with IGB , we've included scripts for making an IGB quickload directory.

     These require some additional sequence files from NCBI. For the example, those can be downloaded by going to the ncbi/ directory and running
     
     sh getSeqs.sh

     Then in the main directory, run

     python3 path-to-xenoGI-directory/misc/createIgbQuckloadDirs.py ncbiIgbDirMap.txt ncbiHumanMap.txt ncbi/ igbExample

     Move species.txt and contents.txt into the directory igbExample.

     Now, running the script

     sh moveInBed.sh
     
     moves the bed files into the newly created igbExample directory. You can then set IGB up to load this.
