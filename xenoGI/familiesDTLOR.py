
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
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation as affinity
from collections import Counter
from .treeParser import *
from .DTLOR_DP import *

import sys
import time

#### Main function

def createDTLORFamiliesO(tree,scoresO,genesO,aabrhHardCoreL,paramD,method="threshold"):
    '''Given a scoresO object, create gene families for the DTLOR approach.
    '''
    synThresholdD = getSynThresholdD(paramD,scoresO,genesO,aabrhHardCoreL,tree)
    initial_families=createInitialFamily(scoresO)
    print("The number of initial families is: ")
    print(len(initial_families))
    locusMap={}
    #construct Family object for the entire problem
    familiesO = Families(tree)
    distribution=[]
    init_dist=[]
    #these locusFamId need to be unique across all
    locusFamId=0
    for index,family in enumerate(initial_families):
        #add each initial family as a Family object (still empty)
        species=list([genesO.numToStrainName(gene) for gene in family])

        mrca=findMRCA(species, tree)
        familiesO.initializeFamily(index,mrca,seedPairL=None) 
        clustered_fam=clusterLocusFamily(family, genesO,scoresO,paramD,synThresholdD,method)
        for locusFamily in clustered_fam: #locusFamily contains all the genes in that lf
            species=[]
            for gene in locusFamily:
                locusMap[str(gene)]=locusFamId
                species.append(genesO.numToStrainName(gene))
            lfMrca=findMRCA(species, tree)
            lf=LocusFamily(index,locusFamId,lfMrca)
            lf.addGenes(locusFamily, genesO)
            familiesO.addLocusFamily(lf)
            locusFamId+=1
        #how many locus families are there in this initial family
        distribution.append(len(clustered_fam))
        init_dist.append(len(family))
    print("The number of locus families is: ")
    print(locusFamId)
       
    distribution=Counter(distribution)
    init_dist=Counter(init_dist)
    # print("Here is the distribution of the number of locus families each initial family has:")
    # print(distribution)
    # print("Distribution of the number of genes in each initial family:")
    # print(init_dist)
    return familiesO, locusMap


    
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

## Supporting functions


def createInitialFamily(scoresO):
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
   
    def stronglyConnected(temp, gene, visited): 
  
       
        visited.add(gene)
        temp.add(gene) 
  
        # Repeat for all vertices adjacent 
        # to this gene
        neighbors=scoresO.getConnectionsGene(gene)
        if neighbors:
            for i in neighbors: 
                if i not in visited:     
                    # Update the list 
                    temp = stronglyConnected(temp, i, visited) 
        
        return temp 

    scoresO.createNodeConnectD() 
    connectedGenes=list(scoresO.nodeConnectD.keys())
    
    initFamily=[]
    visited=set()
    for gene in connectedGenes: 
        if gene not in visited: 
            temp =set()
            newFam=stronglyConnected(temp, gene, visited)
            
            initFamily.append(newFam) 
    
    return initFamily


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

        kmeans=KMeans(n_clusters=num_clusters)
        #cluster the n row vectors into k clusters
        kmeans.fit(V)
        cluster_labels = kmeans.labels_
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
    return a dictionary mapping each gene node to its locus given by the DTLOR recon
    """
    gl_map=dict()
    for key in reconciliation.keys():
        gene=key[0]
        locus_t=key[2]
        locus_b=key[3]
        if gene not in gl_map:
            gl_map[gene]=(locus_t,locus_b)
    return gl_map

def getLeavesInSubtree(startNode, tree_dict):
    if startNode not in tree_dict: #it is a tip
        return [int(startNode)]    #prepare for addGene which needs integer
    else:
        left, right=tree_dict[startNode]
        leaves=getLeavesInSubtree(left, tree_dict)
        leaves.extend(getLeavesInSubtree(right, tree_dict))
        return leaves

def getLocusFamiliesInOrigin(startNode, tree_dict, geneToLocus, origin_num, familiesO, genesO):
    #make sure each node only gets added to the dictionary once
    #get the top and bottom locus of the current gene node branch
    locus_t,locus_b=geneToLocus[startNode]
    if locus_b!=locus_t:
        #new locus family with the current gene node as the MRCA
        lf=LocusFamily(origin_num,int(locus_b),startNode)
        leaves= getLeavesInSubtree(startNode, tree_dict)
        lf.addGenes(leaves, genesO)
        familiesO.addLocusFamily(lf)
    
    if startNode in tree_dict:
        left, right=tree_dict[startNode]
        getLocusFamiliesInOrigin(left, tree_dict, geneToLocus, origin_num, familiesO, genesO)
        getLocusFamiliesInOrigin(right, tree_dict, geneToLocus, origin_num, familiesO, genesO)
        

    
    

def addOriginFamily(reconciliation, geneTree, originFamiliesO, origin_num, genesO):
    geneToLocus=parseReconForLocus(reconciliation)
    tree_dict=getTreeDictionary(geneTree,{})
    #construct an origin family for each Origin event, a locus family for each R event
    for gene in geneToLocus.keys():
        locus_t,locus_b=geneToLocus[gene]
        if locus_t=="*":
            if  locus_b!="*":
                startNode=gene  
                #one gene only visited once for checking O event
                originFamiliesO.initializeFamily(origin_num,startNode,seedPairL=None) 
                #one gene would only be visited once for checking R event as well 
                #since only one ancester can have O event and we are going down from there
                getLocusFamiliesInOrigin(startNode, tree_dict, geneToLocus, origin_num, originFamiliesO, genesO)
                origin_num+=1
    return originFamiliesO, origin_num
    


def reconcile(tree,strainNamesT,scoresO,genesO,aabrhHardCoreL,paramD,method="threshold"):
    #probably read the costs from 
    familiesO, locusMap=createDTLORFamiliesO(tree,scoresO,genesO,aabrhHardCoreL,paramD,method="threshold")
    rawFamilyFN = paramD['rawFamilyFN']
    originFamilyFN=paramD['originFamilyFN']
    geneInfoFN = paramD['geneInfoFN']

    # get the num to name dict, only for strains we're looking at.
    genesO.initializeGeneNumToNameD(geneInfoFN,set(strainNamesT))
    #writing out the initial families to rawFam.out
    f=open(rawFamilyFN,'w')
    for fam in familiesO.iterFamilies():
        f.write(fam.fileStr(genesO)+'\n')
    f.close()
    print("Initial families written out to %s"%rawFamilyFN)
    speciesTree=deepcopy(tree)
    gtFileStem = 'fam'
    workDir = paramD['geneFamilyTreesDir']

    # trees.makeGeneFamilyTrees(paramD,genesO,familiesO) #create a directory and a gene tree for each initial family
    
    allGtFilePath = os.path.join(workDir,gtFileStem+'*.tre')
    D=float(paramD["duplicationCost"])
    T=float(paramD["transferCost"])
    L=float(paramD["lossCost"])
    O=float(paramD["originCost"])
    R=float(paramD["rearrangeCost"])
    print("The parameters used for reconciliation are: D= %.2f, T=%.2f, L=%.2f, O=%.2f, R=%.2f"%(D,T,L,O,R))
    #new families object to record all the origin families
    originFamiliesO = Families(tree)
    origin_num=0
    # for treeFN in list(sorted(glob.glob(allGtFilePath)))[9,14,15]: 
    treeFN= "geneFamilyTrees/fam000004.tre"
    if 1==1:
        #read in unrooted tree
        bpTree = Phylo.read(treeFN, 'newick', rooted=False)
        bpTree.root_at_midpoint()  #root arbitrarily for further rerooting 
        Phylo.draw_ascii(bpTree)
        tabulate_names(bpTree)   #name internal nodes

        #TODO: In the future: Do we need to use tuple format?
        tuple_geneTree=bioPhyloToTupleTree(bpTree)
        phi=getTipMapping(tuple_geneTree,genesO)  #gene to species mapping for these specific gene tree
        all_rootings=get_all_rerootings(tuple_geneTree)
        best_score=float('inf')
        bestMPRs=[]
        start1 = time.time()
        #try all the different rerootings and record the ones and their solutions with the best scores
        print(len(all_rootings))
        for rooting in all_rootings:
            geneTree=rooting
            start2 = time.time()
            MPRs,cost=DP(speciesTree, geneTree, phi, locusMap, D, T, L, O, R)
            end2 = time.time()
            print("The time took to do one reconciliation is: %.4f" %(end2 - start2))
            if cost<best_score: 
                #if the score is better than current best
                #update best score, clear record and add new record
                best_score=cost
                bestMPRs=[]
                bestMPRs.append((geneTree,MPRs))
            elif cost==best_score:
                bestMPRs.append((geneTree,MPRs))
        #sample one MPR from the MPRs for this specific unrooted tree
        end1 = time.time()
        print("The time took to reconcile all the rerootings is: %.4f" %(end1 - start1))
        optGeneRooting,optMPRs=random.choice(bestMPRs)
        optMPR=random.choice(optMPRs)
    
        #update originFamiliesO object
        start = time.time()
        originFamiliesO, origin_num= addOriginFamily(optMPR, geneTree, originFamiliesO, origin_num, genesO)
        end = time.time()
        print("The time took to construct the origins family is: %.2f" %(end - start))
        print("Constructed origin families for %s" %treeFN)

    f=open(originFamilyFN,'w')
    for fam in originFamiliesO.iterFamilies():
        f.write(fam.fileStr(genesO)+'\n')
    f.close()
    print("Origin families written out to %s"%originFamilyFN)
    

    

    # return familiesO, originFamiliesO
    return 


