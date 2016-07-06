"""
  Author: Kevin Heath
  Date: 29 May 2015
  Generated for summer research in the Bush lab """

"""Sample command: python superalign.py -s orthologs.out -n 99 -i nuc.afa -o nuc_aligned.afa 
-p partition.out """


import io, sys, getopt, argparse
from itertools import islice


def main(argv):

    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Input file name (.afa)', required=True)
    requiredNamed.add_argument('-s', '--species', help='orthologs file with species names', required=True)
    requiredNamed.add_argument('-o', '--output', help='Output file name', required=True)
    requiredNamed.add_argument('-n', help='Number of species', type=int, required=True)
    requiredNamed.add_argument('-p', help='Partition file name', required=True)  
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()
    
    global verbose
    verbose = args.verbose

    # Generate a list of species from the ortholog file
    species = getSpecies(args.species, args.n)

    if verbose:
      print 'Species found:'
      for sp in species:
        print sp

    # Create a dictionary with each species and the corresponding genes
    souperAlignment = superAlign(species, args.input)

    # Output the super alignment with all the genes concatenated together as a .fasta file
    writeSuperAlignment(souperAlignment, species, args.output)

    # Output each gene as a partition.
    makePartition(souperAlignment.itervalues().next(), args.p)

    if verbose:
      print "\nWhen cryptography is outlawed, bayl bhgynjf jvyy unir cevinpl."

def getSpecies(filename, numSpecies):
    """Extracts the species names from the orthologs.txt file"""
    if verbose:
      print 'Extracting species...'


    with open(filename, 'r') as f:
      species = []

      while numSpecies > 0:
          line = f.readline()
          temp = line.split()
          species.append(temp[1])
          numSpecies-=1

    return species


def superAlign(species, alignments):
    """Concatenate gene blocks together

    Assumes genes are stored in blocks,
    Each species is listed in the same order"""

    if verbose:
      print 'Generating super alignment...'
    

    numSpecies = len(species)
    
    collection = {sp:[] for sp in species}
    numBlocks = 0
    with open(alignments, 'r') as f:
      while True:
        next_n_lines = list(islice(f, 1, numSpecies*2+1, 2))
        
        if not next_n_lines:
          break
        numBlocks+=1
        for index,key in enumerate(species):
          collection[key].append(next_n_lines[index].strip())

    if verbose:
      print 'Processed '+str(numBlocks)+' blocks'

    return collection



def writeSuperAlignment(superalignment, species, outputfile):
    """Write out super alignment to file"""

    if verbose:
      print 'Writing super alignment, prepare for awesomeness...'

    with open(outputfile, 'w+') as f:
      # for key,value in superalignment.iteritems():
        # f.writelines(['>'+key+'\n', ''.join(value)+'\n'])

      for key in species:
        f.writelines(['>'+key+'\n', ''.join(superalignment.get(key))+'\n'])
      f.write('\n')

def makePartition(geneList, output):
    """Write out the partition.

    Takes a single entry in the superalignment dictionary
    and uses the length of each gene to make the partition/

    Assumes the sequence starts at 1, based on examples"""

    if verbose:
      print 'Writing partition file...'

    start,end = 1,1
    count = 1
    with open(output, 'w+') as f:
      it = iter(geneList)
      gene=''
      while True:
        try:
          gene = it.next()
            
        except StopIteration:
          if verbose:
            print "Exausted list of partitions"
          break
        end=start+len(gene)-1
        temp = "DNA, gene{0}={1}-{2}\n".format(count,start,end)
        f.write(temp)
        start=end+1
        count+=1
        



if __name__ == "__main__":
   main(sys.argv)



    

