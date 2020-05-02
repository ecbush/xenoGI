
#import sys,numpy,random
#from scipy.signal import find_peaks
import os
from .Family import *
from . import scores
from .trees import *
from .genomes import *
from .families import *
from .analysis import printTable
import random
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as LA
import math 
from copy import deepcopy
from scipy.sparse import csgraph
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.cluster import AffinityPropagation as affinity
from collections import Counter
from .treeParser import *
from .DTLOR_DP import *
from collections import deque
import sys
import time

#### Main function

def createDTLORFamiliesO(tree,scoresO,genesO,aabrhHardCoreL,paramD,method="threshold"):
    '''Given a scoresO object, create gene families for the DTLOR approach.
    '''
    #NOTE: temporary threshold for maximum size of gene tree, if above, split with kmeans
    #should note how many intial families were above the threshold initially 
    maxFamilySize=60
    synThresholdD = getSynThresholdD(paramD,scoresO,genesO,aabrhHardCoreL,tree)
    initial_families, degree=createInitialFamily(scoresO, genesO)
    print("Number of initial families before size control")
    print(len(initial_families))
    locusMap={}
    #construct Family object for the entire problem
    familiesO = Families(tree)
    #these locusFamId need to be unique across all, start counting from 1
    
    initial_families.sort(reverse=True,key=len)  #sort by number of genes in descending order
    oversizedFamilies=[family for family in initial_families if len(family)>maxFamilySize]
    print("Number of initial families with more than %d genes: %d"%(maxFamilySize,\
                                len(oversizedFamilies)))
    # for index, family in enumerate(oversized): 
    #     visualize_connectivity(family,"oversized_family_%d"%index,scoresO,degree)
 
    def addFamilyHelper(initial_family,locus_families,familyIndex, locusFamId, genesO, tree,totalAddedTolocusFamilies):
        species=list([genesO.numToStrainName(gene) for gene in initial_family])
        mrca=findMRCA(species, tree)
        #add each initial family as a Family object (still empty)
        familiesO.initializeFamily(familyIndex,mrca) 
        for locusFamily in locus_families: #locusFamily contains all the genes in that lf
            species=[]
            for gene in locusFamily:
                locusMap[str(gene)]=locusFamId
                species.append(genesO.numToStrainName(gene))
            lfMrca=findMRCA(species, tree)
            lf=LocusFamily(familyIndex,locusFamId,lfMrca)
            lf.addGenes(locusFamily, genesO)
            totalAddedTolocusFamilies+=len(locusFamily)
            familiesO.addLocusFamily(lf)
            locusFamId+=1
        familyIndex+=1
        return locusFamId, familyIndex,totalAddedTolocusFamilies

    
    initFam_locusFams=[]
    for i,family in enumerate(oversizedFamilies):
        # visualize_connectivity(family,"oversizedFam_%d"%i,scoresO, degree)
        locus_families=clusterLocusFamily(family, genesO,scoresO,paramD,synThresholdD,method)  #list of lists
        splittings=splitFamByLocus(family, scoresO,maxFamilySize, locus_families)
        for j, (family, locus_families) in enumerate(splittings):
            # visualize_connectivity(family,"%i_postSplit_%d"%(i,j),scoresO)
            initFam_locusFams.append((len(family),family,locus_families))
          
            
    for family in initial_families[len(oversizedFamilies):]:
        locus_families=clusterLocusFamily(family, genesO,scoresO,paramD,synThresholdD,method)
        initFam_locusFams.append((len(family),family,locus_families))

      
    #resort the whole list of initial families by descending size
    #this way we make sure the treeFNs in the directory are also ordered
    initFam_locusFams.sort(reverse=True,key=lambda x: x[0])
    familyIndex=1
    locusFamId=1
    fam_sizes=[]
    totalAddedTolocusFamilies=0
    for size, family, locus_families in initFam_locusFams:
        fam_sizes.append(size)
        locusFamId, familyIndex,totalAddedTolocusFamilies=addFamilyHelper(family,locus_families,familyIndex, locusFamId, genesO, tree,totalAddedTolocusFamilies)
    print("Number of families with 1 gene")
    print(len([x for x in fam_sizes if x==1]))
    print("The number of intial families is: ")
    print(familyIndex)
    print("The number of locus families is: ")
    print(locusFamId)
    print("genes in genesO: %d"%len(set(genesO.iterGenes())))
    print("genes in initial familiesO: %d"%len(familiesO.getAllGenes()))
    return familiesO, locusMap
def splitFamByLocus(family, scoresO,maxFamilySize, locus_families):
    """
    Returns a list of new smaller families with their locus families
    in (family, [locus families]) pair
    """
    origin_size=len(family)
    locus_families.sort(reverse=True,key=len)
    #calculate pairwise RawSc between locus families
    similarity={}
    for i in range(len(locus_families)):
        loc1=locus_families[i]
        if loc1[0] not in similarity:
            similarity[loc1[0]]={}
        for j in range(i+1,len(locus_families)):
            loc2=locus_families[j]
            if loc2[0] not in similarity:
                similarity[loc2[0]]={}
            maxSc=0
            for gene1 in loc1:
                for gene2 in loc2:
                    if scoresO.isEdgePresentByEndNodes(gene1, gene2):
                        rawSc=scoresO.getScoreByEndNodes(gene1, gene2, "rawSc")
                    else: rawSc=0
                    if rawSc>maxSc: maxSc=rawSc
            similarity[loc1[0]][loc2[0]]=maxSc  #use the first gene in the loc as identifier
            similarity[loc2[0]][loc1[0]]=maxSc
    #fill the buckets greedily
    allBuckets=[]
    current_bucket=[]
    currentSize=0
    current_loc=None
    while locus_families!=[]: #while there are still more locus to add
        if current_loc==None:
            current_loc=locus_families[0]  #start with the biggest loc fam in queue
        
        else: #find the loc with the highest similarity score and not in any bucket yet
            #current_loc should already be deleted form the list
            max_ind=np.argmax([similarity[current_loc[0]][other_loc[0]] for other_loc in locus_families])
            current_loc=locus_families[max_ind]
        if len(current_loc)>maxFamilySize:
            print("There is a locus family with size over threshold, putting that as one initial family")
        #try putting into the current bucket
        if (currentSize==0) or (currentSize+len(current_loc)<=maxFamilySize):
            #always allow putting into empty bucket, or if adding to bucket doesn't overflow
            current_bucket.append(current_loc)
            currentSize+=len(current_loc)
            locus_families.remove(current_loc)
        else:#start a new bucket
            allBuckets.append(current_bucket) #archive the current one
            current_bucket=[]
            currentSize=0
            current_loc=None
    allBuckets.append(current_bucket) #archive the current one

    output=[]
    newSize=0
    for bucket in allBuckets: #bucket is a list of locus families (also list)
        newInitialFamily=[gene for locus_family in bucket for gene in locus_family]
        newSize+=len(newInitialFamily)
        output.append((newInitialFamily,bucket))

    return output
def spectral_clustering_family(family,scoresO,maxFamilySize):
    """
    use spectral clustering on the raw BLAST score
     to break a big family further 
     Warning: some cluster might still be above threshold level, need 
     to do this iteratively
    """
    
    num_clusters=math.ceil(len(family)/float(maxFamilySize))
    num_gene=len(family)
    genes=list(family)
    graph=np.identity(num_gene)
    #loop over all pairs of genes
    for i in range(len(genes)):
        for j in range(i+1, len(genes)):
            gene1=genes[i]
            gene2=genes[j]
            if scoresO.isEdgePresentByEndNodes(gene1,gene2):
                graph[i][j]=scoresO.getScoreByEndNodes(gene1, gene2, "rawSc")
                graph[j][i]=graph[i][j]
    sc=SpectralClustering(n_clusters=2,eigen_solver=None,random_state=0,affinity='precomputed',)
    #cluster the n row vectors into k clusters
    sc.fit(graph)
    cluster_labels = sc.labels_
    clustered_fam=[]
    for i in range(num_clusters):
        #get the index of gene in that cluster
        indices=np.argwhere(cluster_labels==i)
        #get the actual index associated with each gene
        fam=set([genes[index[0]] for index in indices])
        clustered_fam.append(fam)
    return clustered_fam

def visualize_connectivity(initial_family,fn,scoresO,degree=None):
    if degree==None: #not the same as the whole gene set 
        degree={}
        genes=list(initial_family)
        for gene in genes:
            all_neighbors=set(scoresO.getConnectionsGene(gene))
            neighbors_in_family=[neighbor for neighbor in genes if neighbor in all_neighbors]
            degree[gene]=len(neighbors_in_family)
    degrees=[degree[gene] for gene in initial_family]
    plt.hist(degrees,bins=20, color='b', edgecolor='k', alpha=0.65)
    plt.axvline(sum(degrees)/float(len(degrees)), color='k', linestyle='dashed', linewidth=1, label='Mean: {:.2f}'.format(sum(degrees)/float(len(degrees))))
    plt.axvline(len(degrees)/float(2), color='c', linestyle='dotted', linewidth=1, label='N/2: {:.2f}'.format(len(degrees)/float(2)))
    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    # min_ylim, max_ylim = plt.ylim()
    # plt.text(sum(degrees)/float(len(degrees))*1.01, max_ylim*0.7, 'Mean: {:.2f}'.format(sum(degrees)/float(len(degrees))))   
    # plt.text(len(degrees)/float(2)*1.01, max_ylim*0.9, 'Half of nodes: {:.2f}'.format(len(degrees)/float(2))) 
    plt.gca().set(title='%s Histogram'%fn, ylabel='frequency', xlabel="Number of neighbors by BLAST")
    if not os.path.exists("initialFamily_histogram"):
        os.mkdir("initialFamily_histogram")

    plt.savefig("initialFamily_histogram/%s"%fn)
    plt.close()


#THIS CURRENTLY RUNS VERY SLOW DUE TO THE STRUCTURE OF THE TREE
def findMRCA(species, tree):
    
    #find the LCA recursively 
    #base cases

    if len(species)==1:
        return species[0]

    specie1=species[0]
    specie2=species[1]

    if len(species)==2:
        return findLCA(specie1,specie2,tree)
    
    else:
        LCA=findLCA(specie1, specie2, tree)

        newSpecies=[LCA]
        newSpecies.extend(species[2:])

        return findMRCA(newSpecies, tree)

  

def findLCA(specie1, specie2, tree):
    if len(tree)==0:
        return None
    elif tree[0]==specie1 or tree[0]==specie2:
        #if root is equal to either, return root as the lCA
        return tree[0]
    else:
        left=findLCA(specie1,specie2, tree[1])
        right=findLCA(specie1, specie2, tree[2])
        if (left is not None) and (right is not None):
            #they are in different subtrees, return the root
            return tree[0]
        if left is None:
            return right
        else: return left




def createInitialFamily(scoresO, genesO):
    '''
    Input
    ------------------------------------------------
    scoresO:        a Score object that scores the 
                    BLAST results from all the genes

    Output
    ------------------------------------------------
    initFamily:      a list of sets where each set stores
                    all the genes that are connected as
                    indicated by significant BLAST score 
                    (that appear in scoresO)
    '''
    degree={}
   
    def stronglyConnected(temp, gene, visited): 
  
       
        visited.add(gene)
        temp.add(gene) 
  
        # Repeat for all vertices adjacent 
        # to this gene
        neighbors=scoresO.getConnectionsGene(gene)
        degree[gene]=len(neighbors)
        if neighbors:
            for i in neighbors: 
                if i not in visited:     
                    # Update the list 
                    temp = stronglyConnected(temp, i, visited) 
        
        return temp 
    allGenes=list(genesO.iterGenes())
    scoresO.createNodeConnectD() 
    connectedGenes=list(scoresO.nodeConnectD.keys())
    initFamily=[]
    visited=set()
    for gene in allGenes: 
        if gene in scoresO.nodeConnectD:
            if gene not in visited: 
                temp =set()
                newFam=stronglyConnected(temp, gene, visited)
                
                initFamily.append(newFam) 
        else:
            fam=set()
            fam.add(gene)
            initFamily.append(fam)
  
    return initFamily, degree


def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)

def consensusEigens(vals1,vecs1, vals2, vecs2):
    vals=np.zeros(vals1.shape)
    vecs=np.zeros(vecs1.shape)
    for i in range(len(vals1)):
        val1=vals1[i]
        val2=vals2[i]
        #ignore the complex part if exist, these are inaccuracy produced by scipy
        if isinstance(val1, complex): val1=val1.real
        if isinstance(val2, complex): val2=val2.real
        #if there is a consensus and the value is greater than 0
       
        greater=max(max(val1,val2),0)
        if greater==val1:
            vals[i]=val1
            vecs[:,i]=vecs1[:,i]
        #if they are all neg, just force to be zero and choose a random eigen vector
        if greater==0:
            vals[i]=0.0
            vecs[:,i]=random.choice([vecs1[:,i],vecs2[:,i]])
        else:
            vals[i]=val2
            vecs[:,i]=vecs2[:,i]
        
    return vals, vecs

#This function is modified from isSameLocusFamily() in families.py
def isSameLocus(gene1,gene2,scoresO,genesO,paramD,synThresholdD):
    
    
    coreSynSc = scoresO.getScoreByEndNodes(gene1,gene2,'coreSynSc')
    synSc = scoresO.getScoreByEndNodes(gene1,gene2,'synSc')
    euclid_score=math.sqrt(coreSynSc**2+synSc**2)/math.sqrt(2)
    if not scoresO.isEdgePresentByEndNodes(gene1,gene2):
        # Within our families, there may be some gene-gene edges
        # missing due to the fact that blast could have just missed
        # significance etc. If the edge isn't there, then we do not
        # have evidence that these genes should be in the same locus
        # family, and we return false.
        
        # NOTE: An alternative approach would be to actually calculate
        # these scores here. But they're likely to be low...

        return (False, euclid_score)
    strain1 = genesO.numToStrainName(gene1)
    strain2 = genesO.numToStrainName(gene2)
    strainPair = tuple(sorted([strain1,strain2]))
    minSynThreshold = synThresholdD['minSynThreshold'][strainPair]
    minCoreSynThreshold = synThresholdD['minCoreSynThreshold'][strainPair]
    
    if coreSynSc < minCoreSynThreshold or synSc < minSynThreshold:
        # one of the two types of synteny below threshold, so this
        # pair doesn't meet the requirements for being in the same
        # LocusFamily
        addIt = False
    else:
        addIt = True
            
    return (addIt, euclid_score)
    
def clusterLocusFamily(family, genesO, scoresO, paramD,synThresholdD,method,min_gene=2, max_cluster=5):
    """
    Input 
    ----------------------------------------
    family:            A initial family of genes generated by 'createInitialFamily'. This is represented 
                            by a set that contains genes that are connected by BLAST results
    scoresO:                The Score object that contains the information for all pairwise synteny scores

    Return 
    ----------------------------------------
    locusFamilies:          A list of lists, sub-lists are each locus family if produced by spectral clustering

    NOTE: Synteny scores are between 0 and 1, giving the percentage of syntenic
            genes shared
    """
    #if there are too few genes, we should probably not cluster further
    if len(family)<=min_gene:
        locusFamily=list(family)
        return [locusFamily]
    #construct the fully connected graph with weighted edges 
    num_gene=len(family)
    genes=np.array(list(family))
    gene_neighbors={}
    #note there is no edge going from one to self.
    if method=="spectral":
        graph=np.zeros((num_gene,num_gene))
    if method=="affinity":
        graph=np.ones((num_gene,num_gene))
    #loop over all pairs of genes
    for i in range(len(genes)):
        for j in range(i+1, len(genes)):
            gene1=genes[i]
            gene2=genes[j]
            if scoresO.isEdgePresentByEndNodes(gene1,gene2):
                sameLocus, euclid_score=isSameLocus(gene1,gene2,scoresO,genesO,paramD,synThresholdD)
                #the graph is symmetrical
                if method=="affinity" or method=="spectral":
                    graph[i][j]=euclid_score
                    graph[j][i]=euclid_score
                if method=="threshold":
                    #initialize the adjacency dictionary if key not already exists
                    if gene1 not in gene_neighbors:
                        gene_neighbors[gene1]=[]
                    if gene2 not in gene_neighbors:
                        gene_neighbors[gene2]=[]
                    #if they are above threshold for both synteny scores, add to each other's neighbors
                    if sameLocus:
                        # graph[i][j]=1
                        # graph[j][i]=1
                        gene_neighbors[gene1].append(gene2)
                        gene_neighbors[gene2].append(gene1)


    
    if method=="spectral":
        #compute the unnormalized Graph Laplacian
        D = np.diag(graph.sum(axis=1))
        laplacian=D-graph
        L=csgraph.laplacian(graph, normed=False)


        vals1, vecs1= LA.eigh(L)
        vals2, vecs2=LA.eig(L)
        vals, vecs=consensusEigens(vals1, vecs1, vals2, vecs2)

        sorted_eigengap_index = np.argsort(np.diff(vals))[::-1]
        #only keep the indices that correspond to cluster less or equal to ax_cluster
        index_largest_gap=list(filter(lambda x: x <= max_cluster-1, sorted_eigengap_index))[0]
        num_clusters = index_largest_gap + 1


        #now, do kmeans
        #let V be the matrix containing the eigenvectors corresponding to k smallest eigenvalues as columns
        smallest_k=np.argsort(vals)[:num_clusters]
        V=vecs[:,smallest_k]

        spectral=KMeans(n_clusters=num_clusters)
        #cluster the n row vectors into k clusters
        spectral.fit(V)
        cluster_labels = spectral.labels_
    if method=="affinity":
        affinity_prop=affinity(convergence_iter=30, affinity="precomputed")
        f = affinity_prop.fit(graph)
        cluster_labels=f.labels_
        if cluster_labels[0]==-1:
            #the algorithm did not converge
            cluster_labels=np.arange(len(genes))
    clustered_fam=[]
    if method!="threshold":
        num_clusters=len(np.unique(cluster_labels))
        for i in range(num_clusters):
            #get the index of gene in that cluster
            indices=np.argwhere(cluster_labels==i)
            #get the actual index associated with each gene
            fam=[genes[index[0]] for index in indices]
            clustered_fam.append(fam)

    if method=="threshold":
        #find connected components
        def stronglyConnected(temp, gene, visited): 
            visited.add(gene)
            temp.append(gene)
            neighbors=gene_neighbors[gene]
            if neighbors:
                for i in neighbors: 
                    if i not in visited:     
                        # Update the list 
                        temp = stronglyConnected(temp, i, visited) 
            return temp 
        #end of helper
        visited=set()
        for gene in genes: 
            if gene not in visited: 
                temp =[]
                newFam=stronglyConnected(temp, gene, visited)
                clustered_fam.append(newFam) 

    

    return clustered_fam


def getTipMapping(geneTree, genesO):
    """
    Fill out the tip mapping (from gene to species) using the binary search function from genomes
    """
    leaves=leafList(geneTree)
    phi={}
    for leaf in leaves:
        geneNum=leaf
        if isinstance(leaf, str):
            geneNum=int(leaf)
        phi[leaf]=genesO.numToStrainName(geneNum)
    return phi


def parseReconForLocus(reconciliation):
    """
    Input
    ---------
    reconciliation:  DTLOR reconciliation dictionary for the whole gene tree
    Return
    ---------
    gl_map:  a dictionary mapping each gene node to its (top, bottom) loci
    given by the DTLOR recon

    mappingNodesWithGene:  a list of keys in the reconciliation that contains the gene mapping
    """
    gl_map=dict()
    mappingNodesWithGene=dict()
    for key in reconciliation.keys():
        gene=key[0]
        locus_t=key[2]
        locus_b=key[3]
        if gene not in gl_map:
            gl_map[gene]=(locus_t,locus_b)
        if gene not in mappingNodesWithGene:
            mappingNodesWithGene[gene]=[key]
        else:
            mappingNodesWithGene[gene].append(key)
    return gl_map, mappingNodesWithGene

def getMRCAforOR(startNode, tree_dict, geneToLocus,MRCA):
    """
    populates the MRCA map from gene leaf to the gene node of the last O or R event
    from it to the root
    """
    locus_t,locus_b=geneToLocus[startNode]
    if locus_b!=locus_t:
        leaves= getLeavesInSubtree(startNode, tree_dict)
        for leaf in leaves:
            MRCA[leaf]=startNode
    if startNode in tree_dict:
        left, right=tree_dict[startNode]
        getMRCAforOR(left, tree_dict, geneToLocus, MRCA)
        getMRCAforOR(right, tree_dict, geneToLocus, MRCA)

    

def getLocusFamiliesInOrigin(startNode, tree_dict, geneToLocus, origin_num, familiesO, genesO, recon):
    """
    Helper function that adds the locus families to originFamiliesO (familiesO). These locus families
    belong in the origin family object that starts with startNode
    """
    LeaftoMRCA={}
    getMRCAforOR(startNode,tree_dict,geneToLocus,LeaftoMRCA)
    MRCAtoLeaves={}
    for gene, MRCA in LeaftoMRCA.items():
        if MRCA in MRCAtoLeaves:
            MRCAtoLeaves[MRCA].append(gene)
        else:
            MRCAtoLeaves[MRCA]=[gene]
    for startNode in MRCAtoLeaves:
        _,locus_b=geneToLocus[startNode]
        _,speciesBr,_=recon[startNode]
        #new locus family with the current gene node's recon mapping species node/branch as the MRCA
        newLocusNum=familiesO.getNumLocusFamilies()
        # lf=LocusFamily(origin_num,int(locus_b),speciesBr)
        lf=LocusFamily(origin_num,newLocusNum,speciesBr, int(locus_b))
        leaves= MRCAtoLeaves[startNode]
        lf.addGenes(leaves, genesO)
        familiesO.addLocusFamily(lf)

def getPartialRecon(reconciliation,mappingNodes):
    """
    Returns a partial dictionary with only key,values pairs
    corresponding to the sub(gene)tree rooted at startNode. 
    Key is a gene node, value is (species branch, locus)
    """
  
    nodeMap={}
    queue=deque(mappingNodes)
    while len(queue)>0:
        currentMapping=queue.popleft()
        value=reconciliation[currentMapping]
        eventType, child1Mapping, child2Mapping=value
        geneBr, speciesBr, tl, bl=currentMapping
        if eventType!='L':  #if it is a loss event the gene node has not been mapped to the bottom 
            nodeMap[geneBr]=(eventType,speciesBr,bl)
        for child in [child1Mapping,child2Mapping]:
            if child!=(None, None, None, None):
                queue.append(child)
    return nodeMap

def addOriginFamily(reconciliation, geneTree,originFamiliesO,genesO, treeFN):
    
    geneToLocus, geneToMappingNodes=parseReconForLocus(reconciliation)
    tree_dict=getTreeDictionary(geneTree,{})
   
    #construct an origin family for each Origin event, a locus family for each R event
    for gene in geneToLocus.keys():
        locus_t,locus_b=geneToLocus[gene]
        if locus_t=='*':
            if  locus_b!='*':  #if you use integer conversion, '*' need to be changed
                startNode=gene  
                #parse the reconciliation to only include the part responsible for origin family
                recon=getPartialRecon(reconciliation, geneToMappingNodes[startNode])
                origin_num=originFamiliesO.getNumFamilies()
                _,speciesBr,_=recon[startNode]
                originFamiliesO.initializeFamily(origin_num,speciesBr,geneTree,recon) #use the species for MRCA
                getLocusFamiliesInOrigin(startNode, tree_dict, geneToLocus, origin_num, originFamiliesO, genesO, recon)
   
 
def runDTLORGroup(argT):
    """
    Run reconciliation a group of tree files
    """
    treeFN_list,speciesTree,locusMap,genesO, D, T, L, O, R= argT
    outputL = []
    for treeFN in treeFN_list:
        optGeneRooting,optMPR = runDTLOR(treeFN,speciesTree,locusMap,genesO, D, T, L, O, R )
        outputL.append((optGeneRooting,optMPR, treeFN))
    return outputL

def runDTLOR(treeFN,speciesTree, locusMap,genesO, D, T, L, O, R ):
    bpTree = Phylo.read(treeFN, 'newick', rooted=False)
    try:
        bpTree.root_at_midpoint()  #root arbitrarily for further rerooting 
    except:
        bpTree.root_with_outgroup(bpTree.get_terminals()[0])
    tabulate_names(bpTree)   #name internal nodes
    if is_binary(bpTree):
        pass
    else: print(treeFN)
    assert is_binary(bpTree)
    locus_map=rerootingPruning(bpTree, locusMap)
    # print("gene to locus map")
    # print([(term.name, locus_map[term.name]) for term in bpTree.get_terminals()])

    # for term in bpTree.get_terminals(): term.name=term.name+"_"+str(locus_map[term.name])  #rename the tips adding the syntenic loc
    # Phylo.draw_ascii(bpTree)
    tuple_geneTree=bioPhyloToTupleTree(bpTree)
    phi=getTipMapping(tuple_geneTree,genesO)  #gene to species mapping for these specific gene tree
    
    all_rootings=get_all_rerootings(tuple_geneTree, locus_map)
    if all_rootings==[]:  #all rerooting not valid (all nodes have the same loc)
        all_rootings=[tuple_geneTree]
    best_score=float('inf')
    bestMPRs=[]
    start1 = time.time()
    #try all the different rerootings and record the ones and their solutions with the best scores
    # print("The number of rerootings is %d" %len(all_rootings))
    for rooting in all_rootings:
        geneTree=rooting
        # start2 = time.time()
        MPR,cost=DP(speciesTree, geneTree, phi, locusMap, D, T, L, O, R)
        # end2 = time.time()
        # print("The time took to do one reconciliation is: %.4f" %(end2 - start2))
        if cost<best_score: 
            #if the score is better than current best
            #update best score, clear record and add new record
            best_score=cost
            bestMPRs=[]
            bestMPRs.append((geneTree,MPR))
        elif cost==best_score:
            bestMPRs.append((geneTree,MPR))
    #sample one MPR from the MPRs for this specific unrooted tree
    end1 = time.time()
    # print("The time took to reconcile all the rerootings is: %.4f" %(end1 - start1))
    optGeneRooting,optMPR=random.choice(bestMPRs) 

    return optGeneRooting,optMPR
    

def reconcile(tree,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,method="threshold"):
    #probably read the costs from 
    
    familiesO, locusMap=createDTLORFamiliesO(tree,scoresO,genesO,aabrhHardCoreL,paramD,method="threshold")
    rawFamilyFN = paramD['rawFamilyFN']
    originFamilyFN=paramD['originFamilyFN']
    geneInfoFN = paramD['geneInfoFN']
    numProcesses = paramD['numProcesses']

    # get the num to name dict, only for strains we're looking at.
    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT))
    #writing out the initial families to rawFam.out
    writeFamilies(familiesO,rawFamilyFN,genesO,strainNamesT,paramD)
    print("Initial families written out to %s"%rawFamilyFN)
    speciesTree=deepcopy(tree)
    gtFileStem = 'fam'
    workDir = paramD['geneFamilyTreesDir']

    trees.makeGeneFamilyTrees(paramD,genesO,familiesO) #create a directory and a gene tree for each initial family 
    
    print("Finished making gene trees")
    allGtFilePath = os.path.join(workDir,gtFileStem+'*.tre')
    D=float(paramD["duplicationCost"])
    T=float(paramD["transferCost"])
    L=float(paramD["lossCost"])
    O=float(paramD["originCost"])
    R=float(paramD["rearrangeCost"])
    print("The parameters used for reconciliation are: D= %.2f, T=%.2f, L=%.2f, O=%.2f, R=%.2f"%(D,T,L,O,R))
    #new families object to record all the origin families
    originFamiliesO = Families(tree)
    #add single-gene initial family as their own origin family
    for init_family in familiesO.iterFamilies():
        familyT = tuple(init_family.iterGenes())  
        if len(familyT) == 1:
            newFamNum=originFamiliesO.getNumFamilies()
            originFamiliesO.initializeFamily(newFamNum,init_family.mrca)
            for locusFamily in init_family.getLocusFamilies():
                locusFamily.famNum=newFamNum
                locusFamily.locusFamNum=originFamiliesO.getNumLocusFamilies()
                originFamiliesO.addLocusFamily(locusFamily)
         
    allTreeFN=list(sorted(glob.glob(allGtFilePath)))
    binaryTreeFN=[]
    no=0
    yes=0
    for treeFN in allTreeFN:
        bpTree = Phylo.read(treeFN, 'newick', rooted=False)
        try:
            bpTree.root_at_midpoint()  #root arbitrarily for further rerooting 
        except:
            bpTree.root_with_outgroup(bpTree.get_terminals()[0])
        tabulate_names(bpTree)   #name internal nodes
        if not is_binary(bpTree): 
            #get fam num
            prefix=workDir+"/"+gtFileStem
            suffix=".tre"
            famNum=treeFN.replace(prefix,'')
            famNum=famNum.replace(suffix,'')
            famNum=int(famNum.lstrip('0'))
            init_fam=familiesO.getFamily(famNum)
            newFamNum=len(originFamiliesO.familiesD)
            tuple_geneTree=bioPhyloToTupleTree(bpTree)
            originFamiliesO.initializeFamily(newFamNum,init_fam.mrca,tuple_geneTree)
            for locusFamily in init_fam.getLocusFamilies():
                locusFamily.famNum=newFamNum
                locusFamily.locusFamNum=originFamiliesO.getNumLocusFamilies()
                originFamiliesO.addLocusFamily(locusFamily)
            no+=1
        else:
            yes+=1
            binaryTreeFN.append(treeFN)
    print("Number of binary trees")
    print(len(binaryTreeFN))
    print("Done adding origin families from nonbinary trees")

    # make list of sets of arguments to be passed to p.map. There
    # should be numProcesses sets.
    argumentL = [([],speciesTree,locusMap,genesO, D, T, L, O, R) for i in range(numProcesses)]
    #distribute all gene tree files into numProcesses separate processes
    for i,treeFN in enumerate(binaryTreeFN):
        argumentL[i%numProcesses][0].append(treeFN)
    
    recon_out=[]
    with Pool(processes=numProcesses) as p:
        # store the results to scoresO as they come in
        for outputL in p.imap_unordered(runDTLORGroup, argumentL):  #each process generates a output list for a group of trees
            recon_out.append(outputL)

    print("finished multiprocessing")
    recon_out=recon_out[::-1]  #reversing the order changes which trees are not adding, which means multiprocess is not where things are wrong
    for outputL in recon_out:
        for optGeneRooting, optMPR, treeFN in outputL:
            #update originFamiliesO object
            addOriginFamily(optMPR, optGeneRooting, originFamiliesO, genesO, treeFN) #we know this is called

    print("Number of genes in origin familiesO: %d"%len(originFamiliesO.getAllGenes()))
    writeFamilies(originFamiliesO,originFamilyFN,genesO,strainNamesT,paramD)
    print("Origin families written out to %s"%(originFamilyFN))

    return 


