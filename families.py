# Functions for a modiefied version of the PhiGs algorithm  http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1523372/
# (We've added the use of synteny)
# Eliot Bush 7/2016
import sys,networkx
import trees

def createSimilarityGraph(scoresFN,geneName2NumD):
    '''Read distances from the distances file and use to create network
with genes and nodes and edges representing global alignment score
between proteins with significant similarity.'''

    G=networkx.Graph()
    for i in range(len(geneName2NumD)): G.add_node(i)
    
    f = open(scoresFN,'r')

    while True:
        s = f.readline()
        if s == '':
            break
        g1,g2,sc=s.split('\t')
        sc = float(sc)

        G.add_edge(geneName2NumD[g1],geneName2NumD[g2],score=sc)

    return G

def nodeDetails(tree,timeOnBrLeadingHere):
    '''Given a tree, return list of (time,node #, left node #, right node
#) for each node. timeOnBrLeadingHere gives the time which should be
added to the time values in the tuple we're making.
    '''
    if tree[1] == ():
        return [(tree[3] + timeOnBrLeadingHere , tree[0], None, None)]
    else:
        l = nodeDetails(tree[1], tree[3] + timeOnBrLeadingHere)
        r = nodeDetails(tree[2], tree[3] + timeOnBrLeadingHere)
        return l + r + [(tree[3] + timeOnBrLeadingHere , tree[0], tree[1][0], tree[2][0])]

        
def createNodeProcessOrderList(tree):
    '''Given a tree, output a list specifying the order of nodes we should
process in the PhiGs algorithm. For each node we give a tuple
(time,node #, left node #, right node #)
    '''
    L=nodeDetails(tree,0)
    L.sort()
    return L

def createLRSets(tree,geneNum2NameD,geneName2StrainNumD):
    '''For every gene in our data, put it into one of three sets. Left,
right, or outgroup. Genes in the left set are found in a species on
the left branch of tree.'''

    leftSpeciesS=set(trees.leafList(tree[1]))
    rightSpeciesS=set(trees.leafList(tree[2]))
    
    leftS=set()
    rightS=set()
    outgroupS=set()
    for geneNum in range(len(geneNum2NameD)): # all genes
        strain=geneName2StrainNumD[geneNum2NameD[geneNum]]
        if strain in leftSpeciesS:
            leftS.add(geneNum)
        elif strain in rightSpeciesS:
            rightS.add(geneNum)
        else:
            outgroupS.add(geneNum)

    return(leftS,rightS,outgroupS)

def closestMatch(gene,rightS,simG):
    '''Find the closest match to gene among the genes in the set geneS in the graph simG. '''

    bestGene=None
    bestEdgeScore = -float('inf')
    for edge in simG.edges_iter(gene):
        if edge[1] in rightS:
            if simG.get_edge_data(*edge)['score'] > bestEdgeScore:
                bestEdgeScore =  simG.get_edge_data(*edge)['score']
                bestGene = edge[1]
    return bestEdgeScore, gene, bestGene
    
def createSeedL(leftS,rightS,simG):
    '''Create a list which has the closest match for each gene on the
opposite side of the tree. e.g. if a gene is in tree[1] then we're
looking for the gene in tree[2] with the closest match. For each gene
we get this closest match, put in a list, sort, and return.'''

    seedL=[]

    for gene in leftS:
        seedL.append(closestMatch(gene,rightS,simG))
    for gene in rightS:
        seedL.append(closestMatch(gene,leftS,simG))

    seedL.sort(reverse=True)
    return seedL

def createAdjacencyGraph(simG,geneOrderT):
    '''Create a graph with genes as nodes, and edges representing adjacent genes.'''

    # get same nodes (genes) as sim graph
    adjG=networkx.Graph()
    for node in simG.nodes_iter(): adjG.add_node(node)

    # add adjacencies as edges
    for contigT in geneOrderT:
        if not contigT == None:
            for geneNumT in contigT:
                for i in range(len(geneNumT)-1):
                    adjG.add_edge(geneNumT[i],geneNumT[i+1])

    return adjG

def pairScore(gn1,gn2,simG):
    '''Given a pair of genes, see if there is an edge. If so, return
score, if not, return 0.'''
    data=simG.get_edge_data(gn1,gn2)
    if data == None:
        return 0
    else:
        return data['score']
    
def synScore(gn1,gn2,simG,adjG):
    '''Given two genes, calculate a synteny score for them.
    e.g. if each list has two genes, it will be the max of this
    o  g1 o       o g1 o
    |     |  or      X
    o  g2 o       o g2 o
    where 
    the g's are our starting pair, and the o's are adjacent genes.
    We divide by two to get an average per adjacent gene, and thus
    this score is between 0 and 1.
    '''
    gn1AdjL = [y for x,y in adjG.edges_iter(gn1)]
    gn2AdjL = [y for x,y in adjG.edges_iter(gn2)]
    if gn1AdjL == [] or gn2AdjL == []:
        return 0
    else:
        sc = 0
        if len(gn1AdjL) == 1:
            for iterGn in gn2AdjL:
                sc += pairScore(gn1AdjL[0],iterGn,simG)
        elif len(gn2AdjL) == 1:
            for iterGn in gn1AdjL:
                sc += pairScore(gn2AdjL[0],iterGn,simG)
        else:
            # both gn1AdjL and gn2AdjL have 2
            #print(gn1AdjL,gn2AdjL,gn1,gn2)
            case1 = pairScore(gn1AdjL[0],gn2AdjL[0],simG) + pairScore(gn1AdjL[1],gn2AdjL[1],simG)
            case2 = pairScore(gn1AdjL[0],gn2AdjL[1],simG) + pairScore(gn1AdjL[1],gn2AdjL[0],simG)
            sc = max(case1,case2)
    return sc / 2.0

def createSynScoresGraph(simG,adjG):
    '''Create a graph with genes as nodes, and edges representing the sum
of the scores of syntenic genes (two immediately adjacent genes).'''

    # get same nodes (genes) as sim graph
    synScoresG=networkx.Graph()
    for node in simG.nodes_iter(): adjG.add_node(node)

    # create and edge for every one in simG. Give it weight of sum of
    # scores of adjacent genes
    for gn1,gn2 in simG.edges_iter():
        sc = synScore(gn1,gn2,simG,adjG)
        synScoresG.add_edge(gn1,gn2,score=sc)
    return synScoresG

def getFamily(seed,synThresh,simG,synScoresG):
    '''Depth first search to get a single cluster.'''

    alThresh,g1,g2 = seed
    
    alreadySearchedS = set()
    notYetSearchedS = set([g1,g2])

    while len(notYetSearchedS) > 0:

        newGenesS=set()
        for gene in notYetSearchedS:
            for edge in simG.edges_iter(gene):
                newGene=edge[1]
                if newGene not in alreadySearchedS and newGene not in newGenesS:
                    sc = simG.get_edge_data(*edge)['score']
                    if sc > alThresh:
                        newGenesS.add(newGene)
                    else:
                        synsc = synScoresG.get_edge_data(*edge)['score']
                        if synsc > synThresh:
                            adjsc = sc + (1 - sc) * synsc
                            if adjsc > alThresh:
                                newGenesS.add(newGene)

            alreadySearchedS.add(gene)
        notYetSearchedS = newGenesS
    return alreadySearchedS
                
    
def families(nodeOrderL,subtreeL,geneNum2NameD,geneName2StrainNumD,synThresh,simG,synScoresG):
    '''Given a graph of genes and their similarity scores (simG) find
families using a PhiGs-like algorithm, with synteny also considered.'''

    geneUsedL = [False for x in geneNum2NameD]
    
    familyL=[]
    
    for t,node,lnode,rnode in nodeOrderL:
        if lnode != None:
            # not a tip
            subtree=subtreeL[node]
            leftS,rightS,outgroupS = createLRSets(subtree,geneNum2NameD,geneName2StrainNumD)
            seedL = createSeedL(leftS,rightS,simG)

            for seed in seedL:
                if seed[0] == -float('inf'):
                    # we've gotten to the point in the seed list with
                    # genes having no match on the other branch
                    break
                clust=getFamily(seed,synThresh,simG,synScoresG)
                #print(node,"len clust",len(clust),seed)

                if any((geneUsedL[gene] for gene in clust)):
                    continue
                else:
                    # none of the genes in this cluster used yet
                    for gene in clust:
                        geneUsedL[gene] = True
                    familyL.append((node,clust))

                    #print("cluster added, now have",len(familyL))

    return familyL


def printFamilies(familyL,geneNum2NameD,geneName2StrainNumD,fileName):
    '''Print all gene families, one family per line. We number families in
order in familyL, and then give each gene with no cluster its own
number.
    '''
    genesInMultiGeneFamsS = set()
    f=open(fileName,'w')
    
    famNum = 0
    for node,fam in familyL:
        genesInMultiGeneFamsS.update(fam) # add all genes in fam
        genesStr = "\t".join((geneNum2NameD[gene] for gene in fam))
        print(famNum,node,genesStr,sep='\t',file=f)
        famNum+=1

    multiGeneFamNum=famNum
    print("Number of multigene families",multiGeneFamNum,file=sys.stderr)

    
    for gene in geneNum2NameD: 
        if not gene in genesInMultiGeneFamsS:
            node = geneName2StrainNumD[geneNum2NameD[gene]]
            print(famNum,node,geneNum2NameD[gene],sep='\t',file=f)
            famNum+=1

    print("Number of single gene families",famNum-multiGeneFamNum,file=sys.stderr)
    print("Number of total families",famNum,file=sys.stderr)
    f.close()
    
def writeG(G,geneNum2NameD,fileName):
    '''Given a graph with genes as nodes, write all edges (pairs of genes)
to file in three columns. Gene 1, gene 2 and score.'''

    f=open(fileName,'w')
    
    for gene1Num,gene2Num in G.edges_iter():
        sc = G.get_edge_data(gene1Num,gene2Num)['score']
        print(geneNum2NameD[gene1Num],geneNum2NameD[gene2Num],format(sc,".6f"),sep='\t',file=f)

    f.close()
