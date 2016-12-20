Code for detecting horizontal transfer events in a clade of bacteria.

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

  . networkx
  https://networkx.github.io/
 

## How to use

- You begin with a set of species with known phylogenetic relationships.

- Create a working directory and then:
  . Obtain genbank (gbff) files for your species and put them in a subdirectory.
  . Create a newick format tree representing their relationships. This tree file should follow the format of example.tre. Note that the branch lengths are dummy in the sense that they are not used for anything. What is important is that internal nodes be given names in this tree. (In the example we label them i0, i1 etc.).
  . create a parameter file based on exampleParams.py. This should have accurate relative (or absolute) paths to the genbank files and the tree file.

- You run the code from within your working directory. If your parameter file was called params.py, and if the xenoGI code was in a parallel directory, you would do this:

python3 ../xenoGI/parseGenbank.py params.py
python3 ../xenoGI/runBlast.py params.py
python3 ../xenoGI/xenoGI.py params.py

- parseGenbank.py runs through the genbank files and produces input files that are used by subsequent code.
- runBlast.py does an all vs. all blast of the genes in these strains. The number of processes it will run in parallel is specified by the numThreads parameter in the parameter file.
- xenoGI.py identifies horizontal transfer events. It is also (mostly) parallelized.

(Note that it works best to use bash at the command line for this. We have observed some problems with other shells.)

- subsequent analysis can be done with analysis.py, e.g.:

python3 -i ../xtrans/analysis.py params.py

From within python, you can then run functions such as

printGroupsAtNode('i0')
or
printFam(10,6782)
