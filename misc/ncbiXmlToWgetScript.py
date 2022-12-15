## Takes a list of assemblies from ncbi in xml format, and produces an
## output file with wget commands for download.
## Authors: Jacob Fischer, Tona Gonzalez, Rachael Soh, Eliot Bush
import sys
import xml.etree.ElementTree as ET

## Example

# Say we have run the following query on ncbi's assembly database:
# Vibrionaceae [Organism] AND "Sequence from type" [Filter] AND "complete genome" [assembly level]

# We save the search results to file in xml format
# (e.g. assembly_results.xml). We next must edit assembly_results.xml
# with a text editor. There needs to be a <data> tag at the very
# beginning, and </data> at the very end. (For some reason this seems
# not to be included in what comes down from ncbi). Then we can use this script as follows:

# python3 ncbiXmlToWgetScript.py assembly_results.xml downloadWget.sh GCF

# the final argument specifies whether we'll get GCF (refseq) assemblies or GCA (genbank).

# downloadWget.sh contains wget commands to get all the assemblies. It can be run like this:

# sh downloadWget.sh

## Functions

def createGenomeLinks(xmlFN, outputFN, assemblySource):
    ''' Parses through xml document output from NCBI assembly search
    and creates and stores genome links in outputString '''
    
    tree = ET.parse(xmlFN)

    root = tree.getroot()

    f = open(outputFN,'w')

    if assemblySource == 'GCA':
        number = 45
    else:
        number = 46

    for child in root:
        link = child[number].text
        
        if link == None:
            continue
        
        i = 0
        end = '' 
        copying = False
        while i != len(link):
            if copying:
                end+=link[i]
                i+=1

            else:
                if link[i:i+4] == "GCA_" or link[i:i+4] == "GCF_":
                    end += link[i]
                    copying = True
                    i+=1
                else:
                    i+=1
        link = "wget " + link + '/' + end

        f.write(link)
        f.write('_genomic.gbff.gz\n')
   
                    
    f.close()

## Main
    
if __name__ == "__main__":
    xmlFN = sys.argv[1]
    outputFN = sys.argv[2]
    assemblySource = sys.argv[3]
    assert(assemblySource in ['GCA', 'GCF'])

    createGenomeLinks(xmlFN, outputFN, assemblySource)
