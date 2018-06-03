import sys,os,shutil
sys.path.append(os.path.join(sys.path[0],'..'))
from xenoGI import parameters

# funcs

def loadNcbiIgbDirMap(fn):
    '''Load a file giving mapping between ncbi file name (stems) and the
directory names we'll use in igbQuickLoad. Return as a list of tuples.'''

    L=[]
    f = open(fn,'r')
    while True:
        s=f.readline()
        if s == '':
            break
        if s == '\n':
            # in case of blank lines, e.g. at end
            continue
        stem,igbStem = s.rstrip().split()
        L.append((stem,igbStem))
    return L

def createOneDir(igbMaindir, igbStem, ncbiDir, ncbiStem, humanStem):
    '''Create a subdirectory and put 2bit, genome.txt and ncib .gff annotation in it.'''
    # make igbStem as subdirectory of igbMaindir
    subdir = igbMaindir + '/' + igbStem
    os.mkdir(subdir)

    call = 'faToTwoBit ' + ncbiDir + ncbiStem + '_genomic.fna ' + subdir + '/' + igbStem + '.2bit'
    os.system(call)

    call = 'twoBitInfo ' + subdir + '/' + igbStem + '.2bit ' + subdir + '/genome.txt'
    os.system(call)
    
    # copy over ncbi .gff file
    shutil.copy(ncbiDir+ncbiStem+'_genomic.gff',subdir+'/'+ncbiStem+'_genomic.gff')

    # make annots.xml file
    f=open(subdir+'/annots.xml','w')
    f.write(annotsString(humanStem,ncbiStem))
    f.close()
    
def annotsString(humanStem,ncbiStem):
    '''Return a string for the annots.xml file, with the proper file names
for the annotations inserted.'''

    a="""<files>
      <file name="""
    b="""	title="xenoGI annotations"
            description="Results from xenoGI annotation pipeline"
            load_hint="Whole Sequence"
            label_field="id"
            background="FFFFFF"
            foreground="000000"
            max_depth="10"
            name_size="12"
            show2tracks="false"/>
      <file name="""

    c="""	title="NCBI annotations"
            description="NCBI's annotations"
            load_hint="Whole Sequence"
            label_field="id"
            background="FFFFFF"
            foreground="000000"
            max_depth="10"
            name_size="12"
            show2tracks="false"/>
    </files>"""

    outStr = a+'"'+humanStem+'-island.bed"\n'+b+'"'+ncbiStem+'_genomic.gff"\n'+c
    
    return outStr
    
if __name__ == "__main__":

    ncbiIgbDirMapFN = sys.argv[1]
    ncbiHumanMapFN =  sys.argv[2]
    ncbiDir = sys.argv[3]
    igbMaindir =  sys.argv[4] # the path to the main quickLoadDir we're making

    ncbiHumanMapD = parameters.loadFileNameMapD(ncbiHumanMapFN)
    ncbiIgbDirMapL = loadNcbiIgbDirMap(ncbiIgbDirMapFN)

    os.mkdir(igbMaindir) # make main dir
    
    for ncbiStem,igbStem in ncbiIgbDirMapL:
        humanStem = ncbiHumanMapD[ncbiStem+'_genomic.gbff']
        createOneDir(igbMaindir, igbStem, ncbiDir, ncbiStem, humanStem)
