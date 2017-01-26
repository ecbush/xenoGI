import sys,statistics,os
import trees, genomes, families, groups


def createGroupByStrainD(leafNodesL,strainNum2StrD,groupByNodeL,familyT,geneNames,geneInfoD):
    '''Return a dict keyed by strain name. Values are lists of tuples
    (groupNum, familyL) where familyL is a list of tuples in the group
    present in that strain. Family tuples are (family,[genes in
    family]).
    '''
    groupByStrainD = {}
    for leaf in leafNodesL:
        groupByStrainD[strainNum2StrD[leaf]]=[]

    # loop over every node in tree. examine each group. from each
    # extract the genes present in each strain and put in right entry
    # in groupByStrainD
    for mrcaNum in range(len(groupByNodeL)):
        for gr in groupByNodeL[mrcaNum]:

            # make dict to collect fams for each strain
            tempStrainD = {}
            for leaf in leafNodesL:
                tempStrainD[strainNum2StrD[leaf]]=[]

            for fam in gr.familyL:

                # famGeneT is a tuple, where index is strain
                # number. It only has non zero entries at
                # tips. Located at each index is a tuple of gene (gene
                # count, (tuple of genes)). We access by looping over
                # leaf nodes
                fgT = familyT[fam].famGeneT

                for leaf in leafNodesL:
                    ct,geneT = fgT[leaf]

                    # if there the family has a gene in this strain,
                    # add tuple (family,[genes in family]) to list for
                    # this strain
                    if ct > 0:
                        geneNamesL=[geneNames.numToName(gene) for gene in geneT]
                        tempStrainD[strainNum2StrD[leaf]].append((fam,geneNamesL))


            # now make group tuple (minStart,group, familyL) where
            # minStart is the lowest start coord we've seen in the
            # group, group is the group id, and familyL is the thing
            # in tempStrainD
            for strain in groupByStrainD:
                # only add if the group is present in this strain

                if tempStrainD[strain] != []:
                    #print(gr.id,mrcaNum,strain)
                    chrom,groupMedianMidpoint,groupMin,groupMax = getGroupPositions(tempStrainD[strain],geneInfoD,strainNum2StrD,gr.id,mrcaNum,strain)
                    grT = (chrom,groupMedianMidpoint,groupMin,groupMax,mrcaNum,gr.id,tempStrainD[strain])
                    groupByStrainD[strain].append(grT)
                    
    # sort each list in groupByStrainD by chrom and groupMedianMidpoint
    for strain in groupByStrainD:
        groupByStrainD[strain].sort(key=lambda x: x[:2])
    
    return groupByStrainD

def getGroupPositions(familyL,geneInfoD,strainNum2StrD,grID,mrcaNum,strain):
    '''Given a list of families (from a single group in a single strain),
return its chrom,start,end.
    '''
    chromL=[]
    groupMin=float('inf')
    groupMax=-float('inf')
    geneMidpointL=[]
    for fam,geneL in familyL:
        for gene in geneL:
            commonName,locusTag,descrip,chrom,start,end,strand=geneInfoD[gene]
            chromL.append(chrom)
            start = int(start)
            end = int(end)
            if start<groupMin:
                groupMin=start
            if end>groupMax:
                groupMax=end

            geneMidpointL.append(int((end-start)/2))
                
    # sanity check: all entries in chromL should be same
    if not all((c==chromL[0] for c in chromL)):
        print("Genes in group",grID,"at mrca",strainNum2StrD[mrcaNum],"in strain",strain,"are not all on the same chromosome.",file=sys.stderr)

    groupMedianMidpoint = statistics.median(geneMidpointL)
    
    return chrom,groupMedianMidpoint,groupMin,groupMax
    
def groupToGff(groupT,geneInfoD,tree,strainNum2StrD,scoreForMRCAatRoot,potentialScoresL,counter):
    '''Given a groupT (the values of groupByStrainD are lists of these)
convert into a string suitable for writing in a gff file. Return
this. Note that we're using the score field to color the genes in
IGB. So, we give genes in a group the same score, and give different
scores to adjacent groups. Counter keeps track of how many groups
we've done already.
    '''
    gffL=[]
    chrom,groupMedianMidpoint,groupMin,groupMax,mrcaNum,grNum,familyL = groupT
    groupID = 'group'+str(grNum)

    # create score for coloring groups
    if trees.isRootNode(tree,mrcaNum):
        score = scoreForMRCAatRoot
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
            gffL.append('\t'.join([chrom,'.','gene',start,end,str(score),strand,'.','ID='+gene+';Name='+Name+';gene='+Name+';Note='+groupID+" | mrca_"+strainNum2StrD[mrcaNum] + " | "+descrip]))

    gffStr = '\n'.join(gffL)
    return gffStr
    
def writeStrainGff(groupByStrainD,geneInfoD,tree,strainNum2StrD,strain,gffFileName,scoreForMRCAatRoot,potentialScoresL):
    '''For a given strain, Write gff file.'''
    counter=0
    f=open(gffFileName,'w')
    f.write('##gff-version 3\n')
    for groupT in groupByStrainD[strain]:
        gffStr = groupToGff(groupT,geneInfoD,tree,strainNum2StrD,scoreForMRCAatRoot,potentialScoresL,counter)
        f.write(gffStr+'\n')
        counter+=1
    f.close()

def createAllGffs(groupByStrainD,geneInfoD,tree,strainNum2StrD,gffFilePath,scoreForMRCAatRoot,potentialScoresL):
    
    # create gff directory
    gffDir = params.gffFilePath.split("*")[0]
    os.mkdir(gffDir)

    gffDir = params.gffFilePath.split("*")[0]
    gffExtension = params.gffFilePath.split("*")[1]
    
    for strain in groupByStrainD:
        gffFileName = gffDir+strain+gffExtension
        writeStrainGff(groupByStrainD,geneInfoD,tree,strainNum2StrD,strain,gffFileName,scoreForMRCAatRoot,potentialScoresL)
    
def createPotentialScoresL(mn,mx,stride,offset):
    '''This function creates a list of potential scores. We use score to
color our groups in IGB. We want adjacent groups to look different,
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
    params = __import__(paramFN.replace('.py', ''))


    ## load data structures we'll use below
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(params.treeFN)
    leafNodesL = trees.leafList(tree)
    geneNames = genomes.geneNames(params.geneOrderFN,strainStr2NumD,strainNum2StrD)
    familyT = families.readFamilies(params.familyFN,tree,geneNames,strainStr2NumD)
    groupByNodeL=groups.readGroups(params.groupOutFN,tree,strainStr2NumD)
    geneInfoD = genomes.readGeneInfoD(params.geneInfoFN)    

    
    # get groups organized by strain
    groupByStrainD = createGroupByStrainD(leafNodesL,strainNum2StrD,groupByNodeL,familyT,geneNames,geneInfoD)

    createAllGffs(groupByStrainD,geneInfoD,tree,strainNum2StrD,params.gffFilePath,params.scoreForMRCAatRoot,params.potentialScoresL)
