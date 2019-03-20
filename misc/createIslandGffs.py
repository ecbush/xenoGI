import sys,statistics,os,glob
from urllib import parse
sys.path.insert(0,os.path.join(sys.path[0],'../'))
from xenoGI import trees, genomes, families, islands, parameters, islandBed

def getIslandPositions(familyL,geneInfoD,strainNum2StrD,grID,mrcaNum,strain):
    '''Given a list of families (from a single island in a single strain),
return its chrom,start,end.
    '''
    chromL=[]
    islandMin=float('inf')
    islandMax=-float('inf')
    geneMidpointL=[]
    for fam,geneL in familyL:
        for gene in geneL:
            commonName,locusTag,descrip,chrom,start,end,strand=geneInfoD[gene]
            chromL.append(chrom)
            start = int(start)
            end = int(end)
            if start<islandMin:
                islandMin=start
            if end>islandMax:
                islandMax=end

            geneMidpointL.append(int((end-start)/2))
                
    # sanity check: all entries in chromL should be same
    if not all((c==chromL[0] for c in chromL)):
        print("Genes in island",grID,"at mrca",strainNum2StrD[mrcaNum],"in strain",strain,"are not all on the same chromosome.",file=sys.stderr)

    islandMedianMidpoint = statistics.median(geneMidpointL)
    
    return chrom,islandMedianMidpoint,islandMin,islandMax

def orderedIslandsInStrain(strainName):
    '''returns a list of all the islands in the strain given, in the order they appear'''
    islandsInStrainL = islandByStrainD[strainName]
    chromStartIslandLL = []
    #for each island in the strain, order first by chromosome and then by island start pos
    for islandT in islandsInStrainL:
        chrom,_,start,_,_,islandNum,familyL=islandT
        chromStartIslandLL.append([chrom,start,islandT])
    sortedChromStartIslandLL = sorted(chromStartIslandLL)
    #make an ordered list of just the islands without the other information
    sortedIslandsL = [x[2] for x in sortedChromStartIslandLL]
    return sortedIslandsL


def islandToGff(islandT,geneInfoD,tree,strainNum2StrD,scoreNodeMapD,potentialScoresL,counter):
    '''Given a islandT (the values of islandByStrainD are lists of these)
convert into a string suitable for writing in a gff file. Return
this. Note that we're using the score field to color the genes in
IGB. So, we give genes in a island the same score, and give different
scores to adjacent islands. Counter keeps track of how many islands
we've done already.
    '''
    gffL=[]
    chrom,islandMedianMidpoint,islandMin,islandMax,mrcaNum,grNum,familyL = islandT
    islandID = 'island_'+str(grNum)

    # create score for coloring islands
    if strainNum2StrD[mrcaNum] in scoreNodeMapD:
        score = scoreNodeMapD[strainNum2StrD[mrcaNum]]
    else:
        score = potentialScoresL[counter%len(potentialScoresL)]
    
    # loop over families to get genes
    for fam,geneL in familyL:
        for gene in geneL:
            commonName,locusTag,descrip,chrom,start,end,strand=geneInfoD[gene]
            if commonName != '':
                Name=commonName
            else:
                Name=gene

            escapedDescrip = parse.quote(descrip) # escape some characters
            attributes = 'ID='+gene+';Name='+Name+';gene='+Name+';Note= | '+islandID+" | fam_"+str(fam)+" | mrca_"+strainNum2StrD[mrcaNum] + " | "+escapedDescrip
            
            gffL.append('\t'.join([chrom,'.','gene',start,end,str(score),strand,'.',attributes]))

    gffStr = '\n'.join(gffL)
    return gffStr
    
def writeStrainGff(islandByStrainD,geneInfoD,tree,strainNum2StrD,strain,gffFileName,scoreNodeMapD,potentialScoresL):
    '''For a given strain, Write gff file.'''
    counter=0
    f=open(gffFileName,'w')
    f.write('##gff-version 3\n')
    orderedIslandsInStrainL = orderedIslandsInStrain(strain)
    for islandT in orderedIslandsInStrainL:
        gffStr = islandToGff(islandT,geneInfoD,tree,strainNum2StrD,scoreNodeMapD,potentialScoresL,counter)
        f.write(gffStr+'\n')
        counter+=1
    f.close()

def createAllGffs(islandByStrainD,geneInfoD,tree,strainNum2StrD,gffFilePath,scoreNodeMapD,potentialScoresL):

    # if directory for gffss doesn't exist yet, make it
    gffDir = gffFilePath.split("*")[0]
    if glob.glob(gffDir)==[]:
        os.mkdir(gffDir)

    gffExtension = gffFilePath.split("*")[1]
    
    for strain in islandByStrainD:
        gffFileName = gffDir+strain+gffExtension
        writeStrainGff(islandByStrainD,geneInfoD,tree,strainNum2StrD,strain,gffFileName,scoreNodeMapD,potentialScoresL)
    
def createPotentialScoresL(mn,mx,stride,offset):
    '''This function creates a list of potential scores. We use score to
color our islands in IGB. We want adjacent islands to look different,
and thus get different scores. So we need a list of scores to march
though, where adjacent things are different. This function produces
such a list. Its not intended to be used every time we run, but rather
to make the list, which can then be put inside a parameter
file. Doubtless there's a better way to make this list.
    '''
    L=[]
    for i in range(0,int(stride/offset)):
        st=mn+offset*i
        end=mx
        for j in range(st,end,stride):
            L.append(j)
    return L
    
if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)

    ## Various other parameters (we're not including these in paramD)
    
    # We want to display different islands with different colors. We do
    # this by giving different 'score' values to different islands. Scores
    # can range from 1 to 1000, but we only give certain discrete scores.
    scoreNodeMapD = {'i1':1,'i2':6} # islands with these mrca values always get this score

    # The score for rest is based on the list below. It was made with
    # createIslandGffs.createPotentialScoresL(100,1001,200,50)
    potentialScoresL=[100, 300, 500, 700, 900, 150, 350, 550, 750, 950, 200, 400, 600, 800, 1000, 250, 450, 650, 850]

    gffFilePath = 'gff/*-island.gff'
    
    ## load data structures we'll use below
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])
    leafNodesL = trees.leafList(tree)
    geneNames = genomes.geneNames(paramD['geneOrderFN'],strainStr2NumD,strainNum2StrD)
    familiesO = families.readFamilies(paramD['familyFN'],tree,geneNames,strainStr2NumD)
    islandByNodeL=islands.readIslands(paramD['islandOutFN'],tree,strainStr2NumD)
    geneInfoD = genomes.readGeneInfoD(paramD['geneInfoFN'])    

    ## Run
    
    # get islands organized by strain
    islandByStrainD = islandBed.createIslandByStrainD(leafNodesL,strainNum2StrD,islandByNodeL,familiesO,geneNames,geneInfoD)

    createAllGffs(islandByStrainD,geneInfoD,tree,strainNum2StrD,gffFilePath,scoreNodeMapD,potentialScoresL)
