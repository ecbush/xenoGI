# Eliot Bush and Matt Wilbur
# 7/2016
import sys, networkx, numpy
from tree import *
from genomes import *

## Functions

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

def createAdjGeneList(geneNum,adjG):
    '''Given a numbered gene, and a graph of adjacencies, produce a list
of length 3 with the gene in the middle and its two neighbors. If
there is no neighbor, put a None there.'''
    gene_neighbors = [y for x,y in adjG.edges_iter(geneNum)]

    if len(gene_neighbors) == 2:
        L = [gene_neighbors[0], geneNum, gene_neighbors[1]]
    elif len(gene_neighbors) == 1:
        L = [gene_neighbors[0], geneNum, None]
    else:
        # no neighbors
        L = [None, geneNum, None]
    return L

    
def createSynSimG(simG,adjG,threshold):
    '''Create a modified graph, synSymG, where the scores have been
adjusted according to synteny.  threshold (between 0 and 1), above
which neighboring gene similarities impact the output similarity
score.
    '''

    synSimG=networkx.Graph()
    for node in simG.nodes_iter(): synSimG.add_node(node)

    for gene1Num,gene2Num in simG.edges_iter():

        # Generate 3x3 matrix of similarity scores of (gene1 + neighbors) and (gene2 + neighbors)
        sim_matrix = numpy.zeros((3, 3))

        sim_matrix_rows = createAdjGeneList(gene1Num,adjG)
        sim_matrix_cols = createAdjGeneList(gene2Num,adjG)

        # we want scores on both diagonals
        for i,j in [(0,0), (1,1), (2,2), (0,2), (2,0)]:

            if sim_matrix_rows[i] == None or sim_matrix_cols[i] == None:
                continue
            row_gene = sim_matrix_rows[i]
            col_gene = sim_matrix_cols[j]
            D=simG.get_edge_data(row_gene,col_gene)
            sim_matrix[i,j] = D['score'] if D != None else 0
        
        # Add the new similarity to the new graph
        old_sim = simG.get_edge_data(gene1Num,gene2Num)['score']
        new_sim = compute_new_similarity(old_sim, sim_matrix, threshold)
        synSimG.add_edge(gene1Num,gene2Num,score=new_sim)
    return synSimG

def writeG(synSimG,geneNum2NameD,synScoresFN):
    '''Given a graph with genes as nodes, write all edges (pairs of genes)
to file in three columns. Gene 1, gene 2 and score.'''

    f=open(synScoresFN,'w')
    
    for gene1Num,gene2Num in synSimG.edges_iter():
        sc = synSimG.get_edge_data(gene1Num,gene2Num)['score']
        print(geneNum2NameD[gene1Num],geneNum2NameD[gene2Num],format(sc,".6f"),sep='\t',file=f)

    f.close()
        
## Matt's funcs for modifying similarity scores based on synteny

def compute_new_similarity(old_sim, sim_matrix, threshold):
    """ Method that computes new similarity score between two genes, based on synteny information.

    Input: 
        sim_matrix:
            A 3x3 numpy matrix of similarity scores. Rows correspond to:

                (gene1leftneighbor, gene1, gene1rightneighbor)

            and columns correspond to

                (gene2leftneighbor, gene2, gene2rightneighbor).

            Entries are gene similarity scores betweeen 0 and 1.

        threshold:
            A float between 0 and 1, the above above which neighboring gene similarities impact the
            output similarity score.

    Output: 
        An updated similarity score between gene1 and gene2. Calculated as an interpolation between
        the old score and 1, based on how high above the threshold the similarities of neighboring
        genes are.
    """
    neighb_sims1 = [sim_matrix[0,0], sim_matrix[2,2]]
    neighb_sims2 = [sim_matrix[0,2], sim_matrix[2,0]]
    
    new_sim1 = compute_interp_similarity(old_sim, neighb_sims1, threshold)
    new_sim2 = compute_interp_similarity(old_sim, neighb_sims2, threshold)

    return max(new_sim1, new_sim2, old_sim)


def compute_interp_similarity(old_sim, neighbor_sims, thresh):
    """ Computes an a new similarity interpolated between old_sim and 1."""
    interp_param1 = 0 if neighbor_sims[0] < thresh else (neighbor_sims[0] - thresh) / (1. - thresh)
    interp_param2 = 0 if neighbor_sims[1] < thresh else (neighbor_sims[1] - thresh) / (1. - thresh)
    avg_interp_param = (interp_param1 + interp_param2) / 2.
    new_sim = old_sim * (1 - avg_interp_param) + 1. * avg_interp_param

    return new_sim



## Main

if __name__ == "__main__":

    paramFN=sys.argv[1]
    params = __import__(paramFN.replace('.py', ''))

    tree,strainStr2NumD,strainNum2StrD = readTree(params.treeFN)
    geneName2NumD,geneNum2NameD,geneName2StrainNumD = createGeneDs(params.geneOrderFN,strainStr2NumD)
    subtreeL=createSubtreeL(tree)
    subtreeL.sort()
    geneOrderT=createGeneOrderTs(params.geneOrderFN,geneName2NumD,subtreeL,strainStr2NumD)
    simG = createSimilarityGraph(params.scoresFN,geneName2NumD)

    
    adjG = createAdjacencyGraph(simG,geneOrderT)

    synSimG = createSynSimG(simG,adjG,.8)

    writeG(synSimG,geneNum2NameD,params.synScoresFN)
