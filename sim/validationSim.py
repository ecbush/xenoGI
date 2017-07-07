import sys,os, copy
sys.path.append(os.path.join(sys.path[0],'../'))
import parameters,trees,genomes,families,islands,analysis

def loadLog(logFN,strainStr2NumD):
    '''Load a simulation log, and store the events in a tuple. index is
node. value is all the events that happened on that branch. We assume
the events at each node come in order in the file.'''

    # load it all into a dict keyed by node
    nodeEventD={}
    f=open(logFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.rstrip().split('\t')
        node=L[0]
        evType = L[1]
        # get the genes. if del, hgt, its on set separated by ' '. If dup, it's two such sets.
        geneL=[]
        for geneSetStr in L[2:]:
            geneSetT = tuple((int(gene) for gene in  geneSetStr.split(' ')))
            geneL.append(geneSetT)
        geneT = tuple(geneL)

        if node in nodeEventD:
            nodeEventD[node].append((evType,geneT))
        else:
            nodeEventD[node] = [(evType,geneT)]

    # convert over to a tuple, where index is node number
    logL=[() for i in range(trees.nodeCount(tree))]
    for nodeStr in nodeEventD:
        logL[strainStr2NumD[nodeStr]] = tuple(nodeEventD[nodeStr])

    return tuple(logL)

def addEvent(geneEventD,key,val):
    '''Add key,val to geneEventD. If key not already in, add it, making
value [val]. If it is already in, append val.'''
    if key in geneEventD:
        geneEventD[key].append(val)
    else:
        geneEventD[key] = [val]
    
def createGeneEventD(logT,strainStr2NumD,initialGeneNumber,tree):
    '''Load a simulation log, and store the events dictionary keyed by
(node,gene) tuples. We assume the events at each node come in order in
the file.
    '''

    # put it all into a dict keyed by node
    geneEventD={}

    # the starting core genes are not created in the log file. do that
    # now.
    for gene in range(initialGeneNumber):
        geneEventD[(tree[0],gene)] = [(tree[0],gene)]

    # now go over logT

    for node in range(trees.nodeCount(tree)):
        eventT = logT[node]
        if eventT != None:
            for event in eventT:
                if event[0] == 'del':
                    for delGene in event[1][0]:
                        addEvent(geneEventD,(node,delGene),('del',None))
                if event[0] == 'dup':
                    sourceGeneL,newDupGeneL = event[1]
                    for i,sourceGene in enumerate(sourceGeneL):
                        addEvent(geneEventD,(node,sourceGene),('dup',newDupGeneL[i]))
    return geneEventD

def famGenes(geneNum, subtree, geneEventD, strainNum2StrD):
    '''Returns the strain genes which descend from gene starting from the
root of subtree. Returns list of strings, where each string is
strainStr-geneNum. Note that geneNum can be different then gene here
if a duplication happens. geneEventD is based on the log.'''

    # process for branch we are on
    geneL=[geneNum]
    if (subtree[0],geneNum) in geneEventD:
        for event in geneEventD[(subtree[0],geneNum)]:
            if event[0] == 'del':
                geneL.remove(geneNum)
            elif event[0] == 'dup':
                geneL.append(event[1]) # add duplicated gene num
        
    # recursively get tip genes
    if subtree[1] == ():
        outGeneL=[]
        for gene in geneL:
            outGeneL.append(str(strainNum2StrD[subtree[0]])+'-'+str(gene))
        return outGeneL
    else:
        outGeneL=[]
        for gene in geneL:
            outGeneL.extend(famGenes(gene, subtree[1], geneEventD, strainNum2StrD))
            outGeneL.extend(famGenes(gene, subtree[2], geneEventD, strainNum2StrD))
        return outGeneL

def createSimIslandGeneSetByNodeL(logT,tree,geneEventD,strainNum2StrD):
    '''Go thought all hgt events in log and get all genes that survive. Output goes in a list organized by node, where the index is branch number on the tree.'''


    simIslandGeneSetByNodeL = [[] for i in range(trees.nodeCount(tree))]

    for node in range(trees.nodeCount(tree)):
        subtree = trees.subtree(tree,node)
        for event in logT[node]:
            if event[0] == 'hgt':
                islandGeneS = set()
                for geneNum in event[1][0]: # event[1][0] has all genes for this hgt
                    # add genes for this family
                    islandGeneS.update(famGenes(geneNum, subtree, geneEventD, strainNum2StrD))
                simIslandGeneSetByNodeL[node].append(islandGeneS)

    return simIslandGeneSetByNodeL


## xgi side

def createXgiIslandGeneSetByNodeL(tree,geneNames,xgiIslandByNodeL,familyL,minIslandLen):
    '''Given a list with the xgi islands organized by node, return another
list organized by node where the elements are sets, each set
containing all the genes for a given xgi island. We only consider xgi
islands that have at least minIslandLen families in them.
    '''

    xgiIslandGeneSetByNodeL = [[] for i in range(trees.nodeCount(tree))]

    for node in range(trees.nodeCount(tree)):
        for islandO in xgiIslandByNodeL[node]:
            if len(islandO.familyL) >= minIslandLen:
                islandGeneS = createXgiIslandGeneSet(islandO,familyL,geneNames)
                xgiIslandGeneSetByNodeL[node].append(islandGeneS)

    return xgiIslandGeneSetByNodeL

            
def createXgiIslandGeneSet(islandO,familyL,geneNames):
    '''Given an island object, get all the genes it is present in and return as a set.'''

    geneS = set()
    for familyNum in islandO.familyL:
        familyO = familyL[familyNum]
        geneS.update(familyO.getGeneNames(geneNames))
    return geneS

    
def validation(simIslandGeneSetByNodeL,xgiIslandGeneSetByNodeL,strainNum2StrD,focalSubtree):
    '''Calculate FDR and TP rate for xenoGI reconstructing islands. For
each simulation island, pick the single best xenoGI island and use
that.'''

    leafL = trees.leafList(focalSubtree)
    
    totAllSim,tipAllSim,internalAllSim=0,0,0
    totAllXgi,tipAllXgi,internalAllXgi=0,0,0
    totMatch,tipMatch,internalMatch=0,0,0
    totExactOverlapCt,tipExactOverlapCt,internalExactOverlapCt=0,0,0
    totSimIslands,tipSimIslands,internalSimIslands=0,0,0
    rowL=[]
    rowL.append(['Node','True pos rate','Pos pred val', 'In single island'])
    rowL.append(['----','-------------','------------', '----------------'])
    tempRowL=[]
    for node in range(trees.nodeCount(focalSubtree)):
        numAllSim,numAllXgi,numMatch,exactOverlapCt=countSimXgiGenesOverlaps(simIslandGeneSetByNodeL[node],xgiIslandGeneSetByNodeL[node])

        totAllSim+=numAllSim
        totAllXgi+=numAllXgi
        totMatch+=numMatch
        totExactOverlapCt+=exactOverlapCt
        totSimIslands += len(simIslandGeneSetByNodeL[node])

        if node in leafL:
            tipAllSim+=numAllSim
            tipAllXgi+=numAllXgi
            tipMatch+=numMatch
            tipExactOverlapCt+=exactOverlapCt
            tipSimIslands += len(simIslandGeneSetByNodeL[node])
        else:
            internalAllSim+=numAllSim
            internalAllXgi+=numAllXgi
            internalMatch+=numMatch
            internalExactOverlapCt+=exactOverlapCt
            internalSimIslands += len(simIslandGeneSetByNodeL[node])
            
        addRow(tempRowL,strainNum2StrD[node],numAllSim,numAllXgi,numMatch,exactOverlapCt,len(simIslandGeneSetByNodeL[node]))

    addRow(rowL,"All nodes",totAllSim,totAllXgi,totMatch,totExactOverlapCt,totSimIslands)
    addRow(rowL,"Tip nodes",tipAllSim,tipAllXgi,tipMatch,tipExactOverlapCt,tipSimIslands)
    addRow(rowL,"Int nodes",internalAllSim,internalAllXgi,internalMatch,internalExactOverlapCt,internalSimIslands)
    rowL.append(['    ','             ','            ','                    '])
    rowL.extend(tempRowL)
    analysis.printTable(rowL,indent=2,fileF=sys.stdout)

def countSimXgiGenesOverlaps(simSL,xgiSL):
    '''Given a list of gene sets corresponding to the simulation islands
(at a node) and a similar list from xenoGI, count the total number of
sim genes, xgi genes, and overlapping genes. Returns total num sim
genes, total num xgi genes, total num overlap genes (this allows more
than one xgi island to be counted), and the count of sim islands with
exactly one xgi island matching.

    '''


    # get all genes of each type
    allSimGenesS=set()
    for S in simSL:
        allSimGenesS.update(S)

    allXgiGenesS=set()
    for S in xgiSL:
        allXgiGenesS.update(S)

    matchS = allSimGenesS.intersection(allXgiGenesS)
        
    # get prop sim islands with exactly 1 overlap among xgi islands
    exactlyOneOverlapCt = 0
    for simS in simSL:
        numOverlap=0
        for xgiS in xgiSL:
            if len(simS.intersection(xgiS)) > 0:
                numOverlap +=1
        if numOverlap == 1:
            exactlyOneOverlapCt += 1
            
    return len(allSimGenesS),len(allXgiGenesS),len(matchS),exactlyOneOverlapCt
    
def addRow(rowL,nodeStr,allSim,allXgi,match,exactOverlapCt,numSimIslands):
    '''Add a row to rowL with true positive rate and positive predictive
value. Returns None.'''

    tpr = format(match/allSim,".3f") if allSim>0 else "None"
    ppv = format(match/allXgi,".3f") if allXgi>0 else "None"
    prop = format(exactOverlapCt / numSimIslands,".3f") if numSimIslands>0 else "None"
    rowL.append([nodeStr,tpr,ppv,prop])


## main

if __name__ == "__main__":

    simParamFN=sys.argv[1]
    simParamD = parameters.loadParametersD(simParamFN)

    xgiParamFN=sys.argv[2]
    xgiParamD = parameters.loadParametersD(xgiParamFN)

    minIslandLen = 2 # min number of families that must be in an xgi
                     # island to count
    
    tree,strainStr2NumD,strainNum2StrD = trees.readTree(simParamD['treeFN'])
    geneNames = genomes.geneNames(xgiParamD['geneOrderFN'],strainStr2NumD,strainNum2StrD)

    
    logT = loadLog(simParamD['logFile'],strainStr2NumD)

    geneEventD=createGeneEventD(logT,strainStr2NumD,simParamD['initialGeneNumber'],tree)

    # get sets of genes from each island
    simIslandGeneSetByNodeL = createSimIslandGeneSetByNodeL(logT,tree,geneEventD,strainNum2StrD)

    

    # load xgi families,islands
    familyL = families.readFamilies(xgiParamD['familyFN'],tree,geneNames,strainStr2NumD)
    xgiIslandByNodeL=islands.readIslands(xgiParamD['islandOutFN'],tree,strainStr2NumD)

    # get sets of genes for each island
    xgiIslandGeneSetByNodeL = createXgiIslandGeneSetByNodeL(tree,geneNames,xgiIslandByNodeL,familyL,minIslandLen)

    focalSubtree=trees.subtree(tree,strainStr2NumD[xgiParamD['rootFocalClade']])
    validation(simIslandGeneSetByNodeL,xgiIslandGeneSetByNodeL,strainNum2StrD,focalSubtree)
