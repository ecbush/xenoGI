from collections import defaultdict
from enum import Enum, auto
from itertools import product
from functools import reduce
from .DTLOR_DP import check_tip, preorder, postorder, delta_r, find_min_events

Infinity = float("inf")

class GraphType(Enum):
    """
    Defines how each NodeType interacts with the graph when
    choosing a reconciliation. Choose is for "mapping" nodes, that
    choose which optimal event below them to use. All is for "event nodes",
    all of whose children must be used if the event is chosen.
    """
    CHOOSE = auto()
    ALL = auto()

class NodeType(Enum):
    """
    Encodes the type of a node in a reconciliation graph
    """
    # DTL events
    COSPECIATION = (1, GraphType.ALL)
    DUPLICATION = (2, GraphType.ALL)
    LOSS = (3, GraphType.ALL)
    TRANSFER = (4, GraphType.ALL)
    # Assigns locations to the left and right children of a gene node
    LOCATION_ASSIGNMENT = (5, GraphType.ALL)
    # Maps a gene tree node to a list of locations
    LOCATION_LIST = (6, GraphType.CHOOSE)
    # Maps a gene tree node to a syntenic location
    LOCATION_MAPPING = (7, GraphType.CHOOSE)
    # Maps a gene tree node to a list of species tree nodes
    # Used for choosing the optimal root of a reconciliation
    SPECIES_LIST = (8, GraphType.CHOOSE)
    # Maps gene tree node to species tree node
    SPECIES_MAPPING = (9, GraphType.CHOOSE)
    # When a gene gets a true location (which location is yet unknown)
    ORIGIN = (10, GraphType.ALL)
    # "Root" event for choosing the syntenic location of the root
    ROOT = (11, GraphType.CHOOSE)
    # Event node that represents a change from one syntenic location to another
    REARRANGEMENT = (12, GraphType.ALL)
    # Event node that represents the origin of an actual syntenic location
    ORIGIN_EVENT = (13, GraphType.ALL)
    def __init__(self, e_id, graph_type):
        self._e_id = e_id
        self.graph_type = graph_type
    def __repr__(self):
        return "{}".format(self.name)
    def __eq__(self, other):
        return self._e_id == other._e_id
    def __ne__(self, other):
        return self._e_id != other._e_id
    def __hash__(self):
        return hash((self._e_id))

# Weights for each event type for the distance metric
default_event_weights = {t: 1 for t in NodeType}
# Set of nodes which are implicated in the distance metric
# For the location mapping node symmetric set distance
dist_matters_nodes = set([
    NodeType.LOCATION_MAPPING,
    NodeType.DUPLICATION,
    NodeType.LOSS,
    NodeType.TRANSFER,
    NodeType.COSPECIATION])
# For the event symmetric set distance
dist_matters_events = set([
    NodeType.ORIGIN_EVENT,
    NodeType.REARRANGEMENT,
    NodeType.DUPLICATION,
    NodeType.LOSS,
    NodeType.TRANSFER,
    NodeType.COSPECIATION])

def compute_dtlor_graph(species_tree, gene_tree, phi, locus_map, D, T, L, O, R):
    gene_root = next(iter(gene_tree))
    species_root = next(iter(species_tree))
    # First, compute C
    C, C_star, C_graph = DTL_reconcile(species_tree, gene_tree, phi, D, T, L)
    S, S_star, S_graph = synteny_reconcile(species_tree, gene_tree, locus_map, R)
    #print(gene_tree)
    #print(S_graph)
    #print(postorder(gene_tree))
    # Union the two graphs before adding Null events
    G = {**C_graph, **S_graph}
    Origin = {}
    Null = {}
    # Compute the Origin and Null tables
    for eg in postorder(gene_tree):
        _, vp, eg1, eg2 = gene_tree[eg]
        o_cost = C_star[eg][0] + S_star[eg][0] + O
        # Each origin event for gene g
        # leads to a species list node for optimal places to put that gene
        c_choice = (NodeType.SPECIES_LIST, eg)
        # And a location list node for optimal syntenic locations to give it
        s_choice = (NodeType.LOCATION_LIST, eg)
        o_event = (NodeType.ORIGIN, eg)
        G[o_event] = [c_choice, s_choice]
        Origin[eg] = o_cost, [o_event]
        if check_tip(vp, eg1, eg2):
            Null[eg] = Infinity
        else:
            # Null or Origin for left child
            left_null = Null[eg1], [(NodeType.LOCATION_MAPPING, eg1, "*")]
            left_origin = Origin[eg1]
            left_cost, left_nodes = find_min_events([left_null, left_origin])
            # Null or Origin for right child
            right_null = Null[eg2], [(NodeType.LOCATION_MAPPING, eg2, "*")]
            right_origin = Origin[eg2]
            right_cost, right_nodes = find_min_events([right_null, right_origin])
            a_nodes = []
            # Location assignment nodes assign both left and right
            for l_node, r_node in product(left_nodes, right_nodes):
                node = (NodeType.LOCATION_ASSIGNMENT, l_node, r_node)
                G[node] = [l_node, r_node]
                a_nodes.append(node)
            cost = left_cost + right_cost
            Null[eg] = cost
            G[(NodeType.LOCATION_MAPPING, eg, "*")] = a_nodes
    # Compute the final choice node and cost: does the root get a location or "*"?
    root_null = (Null[gene_root], [(NodeType.LOCATION_MAPPING, gene_root, "*")])
    root_origin = Origin[gene_root]
    root_cost, root_events = find_min_events([root_null, root_origin])
    G[(NodeType.ROOT,)] = root_events
    #print(G)
    # Remove the non-optimal parts of G
    G = prune_graph(G)
    return root_cost, G

def DTL_reconcile(species_tree, gene_tree, phi, D, T, L):
    A = {}
    # node -> cost
    C = {}
    # node -> [node]
    C_graph = {}
    O = {}
    # Holds the optimal species nodes to map a given gene node
    # gene_node -> (cost, [species_node])
    C_star = {}
    best_switch = {}
    species_root = next(iter(species_tree))
    for eg in postorder(gene_tree):
        C_star[eg] = (Infinity, [])
        _, vp, eg1, eg2 = gene_tree[eg]
        vp_is_a_tip = check_tip(vp, eg1, eg2)
        for es in postorder(species_tree):
            _, vh, es1, es2 = species_tree[es]
            vh_is_a_tip = check_tip(vh, es1, es2)
            # Relevant mapping nodes, for convenience
            eg_es_m = (NodeType.SPECIES_MAPPING, eg, es)
            eg_es1_m = (NodeType.SPECIES_MAPPING, eg, es1)
            eg_es2_m = (NodeType.SPECIES_MAPPING, eg, es2)
            eg1_es_m = (NodeType.SPECIES_MAPPING, eg1, es)
            eg2_es_m = (NodeType.SPECIES_MAPPING, eg2, es)
            # Compute A
            if vh_is_a_tip:
                # Must match the tip mapping phi
                if vp_is_a_tip and phi[eg] == es:
                    A[eg_es_m] = (0, [])
                else:
                    A[eg_es_m] = (Infinity, [])
            else:
                # Cospeciation
                if not vp_is_a_tip:
                    co1_event = (NodeType.COSPECIATION,
                            (NodeType.SPECIES_MAPPING, eg1, es1),
                            (NodeType.SPECIES_MAPPING, eg2, es2))
                    co1_cost = C[(NodeType.SPECIES_MAPPING, eg1, es1)] + \
                            C[(NodeType.SPECIES_MAPPING, eg2, es2)]
                    co1 = (co1_cost, [co1_event])
                    co2_event = (NodeType.COSPECIATION,
                            (NodeType.SPECIES_MAPPING, eg1, es2),
                            (NodeType.SPECIES_MAPPING, eg2, es1))
                    co2_cost = C[(NodeType.SPECIES_MAPPING, eg1, es2)] + \
                            C[(NodeType.SPECIES_MAPPING, eg2, es1)]
                    co2 = (co2_cost, [co2_event])
                    cospeciation = find_min_events([co1, co2])
                else:
                    cospeciation = (Infinity, [])
                # Loss
                loss_es1 = (C[eg_es2_m] + L, [(NodeType.LOSS,
                    eg_es2_m, None)])
                loss_es2 = (C[eg_es1_m] + L, [(NodeType.LOSS,
                    eg_es1_m, None)])
                loss = find_min_events([loss_es1, loss_es2])
                A[eg_es_m] = find_min_events([cospeciation, loss])
            # Compute C
            # Duplication
            if not vp_is_a_tip:
                dup_event = (NodeType.DUPLICATION,
                        eg1_es_m, 
                        eg2_es_m)
                duplication = (D + C[eg1_es_m] + C[eg2_es_m], [dup_event])
            else:
                duplication = (Infinity, [])
            # Transfer
            if not vp_is_a_tip:
                eg2_cost, eg2_locations = best_switch[eg2_es_m]
                eg2_switch_cost = T + C[eg1_es_m] + eg2_cost
                # location[2] extracts the species node from the mapping node that stores the location
                eg2_switch_events = [(NodeType.TRANSFER,
                    (NodeType.SPECIES_MAPPING, eg1, es),
                    (NodeType.SPECIES_MAPPING, eg2, location[2])) \
                        for location in eg2_locations]
                eg2_switch = (eg2_switch_cost, eg2_switch_events)
                eg1_cost, eg1_locations = best_switch[eg1_es_m]
                eg1_switch_cost = T + C[eg2_es_m] + eg1_cost
                eg1_switch_events = [(NodeType.TRANSFER, 
                    (NodeType.SPECIES_MAPPING, eg2, es),
                    (NodeType.SPECIES_MAPPING, eg1, location[2])) \
                        for location in eg1_locations]
                eg1_switch = (eg1_switch_cost, eg1_switch_events)
                transfer = find_min_events([eg2_switch, eg1_switch])
            else:
                transfer = (Infinity, [])
            cost, events = find_min_events([A[eg_es_m], duplication, transfer])
            C[eg_es_m] = cost
            C_graph[eg_es_m] = events
            # Nodes the graph have children which are a list
            # Add the optimal events to the graph with their children
            for event in C_graph[eg_es_m]:
                t,l,r = event
                event_nodes = []
                if l is not None:
                    event_nodes.append(l)
                # Can be none for a loss
                if r is not None:
                    event_nodes.append(r)
                C_graph[event] = event_nodes
            C_star[eg] = find_min_events([C_star[eg], (C[eg_es_m], [es])])
            # Compute O: (cost, [mapping_node])
            O_c = (C[eg_es_m], [eg_es_m])
            if vh_is_a_tip:
                O[eg_es_m] = O_c
            else:
                O[eg_es_m] = find_min_events([O_c, O[eg_es1_m], O[eg_es2_m]])
        # Compute best_switch
        best_switch[(NodeType.SPECIES_MAPPING, eg, species_root)] = (Infinity, [])
        for es in preorder(species_tree):
            _, vh, es1, es2 = species_tree[es]
            vh_is_a_tip = check_tip(vh, es1, es2)
            eg_es_m = (NodeType.SPECIES_MAPPING, eg, es)
            eg_es1_m = (NodeType.SPECIES_MAPPING, eg, es1)
            eg_es2_m = (NodeType.SPECIES_MAPPING, eg, es2)
            if not vh_is_a_tip:
                best_switch[eg_es1_m] = find_min_events([best_switch[eg_es_m], O[eg_es2_m]])
                best_switch[eg_es2_m] = find_min_events([best_switch[eg_es_m], O[eg_es1_m]])
        # Now that we're done computing C_star[eg], add appropriate choice events
        species_choice = (NodeType.SPECIES_LIST, eg)
        C_graph[species_choice] = [(NodeType.SPECIES_MAPPING, eg, es) for es in C_star[eg][1]]
    return C, C_star, C_graph

#TODO gene/species vs. gene/species
def synteny_reconcile(species_tree, gene_tree, locus_map, R):
    # Holds the optimal syntenic locations for a given gene node
    # gene_node -> [location]
    S_star = {}
    # (g,s) -> cost
    S = {}
    # node -> [[node]]
    S_graph = {}
    allsynteny = set(locus_map.values())
    # Capture R for convenience
    delta = lambda l1, l2: delta_r(l1, l2, R)
    for eg in postorder(gene_tree):
        S_star[eg] = (Infinity, [])
        _, vp, eg1, eg2 = gene_tree[eg]
        for lp in allsynteny:
            eg_lp_m = (NodeType.LOCATION_MAPPING, eg, lp)
            if check_tip(vp, eg1, eg2):
                if locus_map[eg] == lp:
                    S[eg_lp_m] = 0
                    S_graph[eg_lp_m] = []
                    # lp is always the optimum for eg with no cost
                    S_star[eg] = (0, [lp])
                else:
                    S[eg_lp_m] = Infinity
            else:
                # Synteny cost for the left child
                l_keep_m = (NodeType.LOCATION_MAPPING, eg1, lp)
                l_keep = (S[l_keep_m], [l_keep_m])
                l_rearrange = (S_star[eg1][0] + R, [(NodeType.LOCATION_LIST, eg1)])
                l_cost, l_nodes = find_min_events([l_keep, l_rearrange])
                # Synteny cost for the right child
                r_keep_m = (NodeType.LOCATION_MAPPING, eg2, lp)
                r_keep = (S[r_keep_m], [r_keep_m])
                r_rearrange = (S_star[eg2][0] + R, [(NodeType.LOCATION_LIST, eg2)])
                r_cost, r_nodes = find_min_events([r_keep, r_rearrange])
                a_nodes = []
                # Create the appropriate assignment nodes
                for l_node, r_node in product(l_nodes, r_nodes):
                    node = (NodeType.LOCATION_ASSIGNMENT, l_node, r_node)
                    S_graph[node] = [l_node, r_node]
                    a_nodes.append(node)
                cost = l_cost + r_cost
                #assert cost < Infinity
                S[(NodeType.LOCATION_MAPPING, eg, lp)] = cost
                S_graph[eg_lp_m] = a_nodes
                S_star[eg] = find_min_events([S_star[eg], (cost, [lp])])
        # Create the appropriate choice node
        location_choice = (NodeType.LOCATION_LIST, eg)
        S_graph[location_choice] = [(NodeType.LOCATION_MAPPING, eg, lp) for lp in S_star[eg][1]]
    return S, S_star, S_graph

#TODO: methods of Graph class?
def prune_graph(G):
    """
    Removes the unnecessary nodes from G.
    """
    new_G = {}
    # Use a set to prevent reprocessing
    extant_nodes = set([(NodeType.ROOT,)])
    while len(extant_nodes) != 0:
        node = extant_nodes.pop()
        if node[0].graph_type is not None:
            #print(node)
            children = G[node]
            #print(children)
            new_G[node] = children
            extant_nodes |= set(children)
    return new_G

def find_MPR(G, rand=False):
    """
    Samples a traversal from a graph
    """
    MPR = {}
    # BFS
    extant_nodes = [(NodeType.ROOT,)]
    while len(extant_nodes) != 0:
        node = extant_nodes.pop()
        children = G[node]
        if len(children) > 0:
            # Choose or take all children based on the GraphType
            if node[0].graph_type is GraphType.CHOOSE:
                if rand:
                    choice = random.choice(children)
                else:
                    choice = children[0]
                MPR[node] = [choice]
                extant_nodes.append(choice)
            elif node[0].graph_type is GraphType.ALL:
                MPR[node] = children
                extant_nodes.extend(children)
            else:
                assert False, "Bad GraphType"
    return MPR

def graph_search_order(G):
    """
    Iterator over the nodes of G in BFS order (parents before children)
    Uses set, so the order isn't deterministic, but for the bottom-up/top-down
    DP algorithms, this does not matter
    """
    extant_nodes = set([(NodeType.ROOT,)])
    while len(extant_nodes) != 0:
        node = extant_nodes.pop()
        node_children = G[node]
        extant_nodes |= set(node_children)
        yield node

def count_MPRs(G):
    """
    Count the total number of traversals in G (which is the number of optimal MPRs)
    Computes a node -> count table which is the number of optimal sub-MPRs in
    the subtree rooted at node.
    """
    counts = {}
    # Reverse the graph search order to get children before parents
    postorder = list(graph_search_order(G))[::-1]
    for node in postorder:
        children = G[node]
        # Base case: nodes with no children
        if len(children) == 0:
            counts[node] = 1
        else:
            child_counts = [counts[child] for child in children]
            # All means multiply
            if node[0].graph_type is GraphType.ALL:
                counts[node] = reduce(lambda x, y: x * y, child_counts)
            # Choose means add
            elif node[0].graph_type is GraphType.CHOOSE:
                counts[node] = sum(child_counts)
            else:
                assert False, "Bad GraphType"
    return counts

def event_frequencies(G, counts):
    """
    Compute the 'frequency' of every node - the number of MPRs it appears in
    """
    frequencies = defaultdict(int)
    for node in graph_search_order(G):
        if node[0] is NodeType.ROOT:
            frequencies[node] = counts[node]
        # Now compute the frequencies for this node's children
        children = G[node]
        multiplier = frequencies[node] / counts[node]
        for child in children:
            # All passes the frequency along
            if node[0].graph_type is GraphType.ALL:
                frequencies[child] += frequencies[node]
            # Choose divides the frequency according to the relative proportion of the count
            elif node[0].graph_type is GraphType.CHOOSE:
                frequencies[child] += multiplier * counts[child]
            else:
                assert False, "Bad GraphType"
    return frequencies

def median_subgraph(G, frequencies, dist_matters):
    """
    Compute the medians by maximizing the sum of the (pre-adjusted) frequencies
    """
    median_graph = {}
    freq_sums = {}
    # Children to parents
    postorder = list(graph_search_order(G))[::-1]
    for node in postorder:
        children = G[node]
        if len(children) == 0:
            freq_sums[node] = 0
            median_graph[node] = []
        else:
            # Choose the score from the best child
            if node[0].graph_type is GraphType.CHOOSE:
                # Negate them because we want to maximize the freq sum
                child_freq_sums = [(-freq_sums[c], [c]) for c in children]
                # Finding the min events of the negative freq_sums finds the maximum freq_sums
                neg_freq_sum, opt_children = find_min_events(child_freq_sums)
                freq_sums[node] = -neg_freq_sum
                median_graph[node] = opt_children
            # Add the freq_sums of both children and itself
            elif node[0].graph_type is GraphType.ALL:
                median_graph[node] = children
                freq_sums[node] = sum([freq_sums[c] for c in children])
            else:
                assert False, "Bad GraphType"
        if node[0] in dist_matters:
            freq_sums[node] += frequencies[node]
    # Need to prune it because it was built bottom-up, so many suboptimal parts have been added
    median_graph = prune_graph(median_graph)
    return median_graph

#TODO: event weights
def build_median_graph(G, event_weights, dist_matters):
    """
    Build the median graph by doing the 3 necessary DP algorithms
    G: Reconciliation graph
    event_weights: dict from NodeType -> float indicating how much to weight each event
    dist_matters: set of NodeType which indicates which types of nodes are included in the distance metric
    """
    # First, get the counts and node frequencies
    counts = count_MPRs(G)
    freqs = event_frequencies(G, counts)
    # Adjust the frequencies by half to get a median
    adj = 0.5 * counts[(NodeType.ROOT,)]
    adjusted_freqs = {node: event_weights[node[0]] * (freq - adj) for node,freq in freqs.items()}
    return median_subgraph(G, adjusted_freqs, dist_matters)

def build_node_median_graph(G, event_weights=default_event_weights):
    """
    Build the median graph for the location mapping node symmetric set distance
    """
    build_median_graph(G, event_weights, dist_matters_nodes)

def build_event_median_graph(G, event_weights=default_event_weights):
    """
    Build the median graph for the event symmetric set distance
    """
    event_graph = build_event_graph(G)
    return build_median_graph(event_graph, event_weights, dist_matters_events)

def get_mapping_nodes(location_nodes, G):
    """
    Get the list of location mapping nodes below  the given location assignments
    """
    maps = []
    for l in location_nodes:
        if l[0] is NodeType.LOCATION_MAPPING:
            maps.append(l)
        elif l[0] is NodeType.LOCATION_LIST:
            maps.extend(G[l])
        else:
            assert False, "Bad node {}".format(l)
    return maps

def create_r_events(node, G, event_graph):
    """
    Take a location mapping node and create the set of R event nodes and location
    assignments that should appear below it (that replace the location list nodes)
    """
    assignments = []
    children = G[node]
    # Create an assignment node for the left and for the right
    for i in [1,2]:
        nodes = set([location_assignment[i] for location_assignment in children])
        maps = get_mapping_nodes(nodes, G)
        if len(maps) > 0:
            assignment = (NodeType.LOCATION_LIST, node, i-1)
            assignments.append(assignment)
            event_graph[assignment] = []
            for m in maps:
                # If the syntenic location matches, it's not an R event
                if m[2] == node[2]:
                    event_graph[assignment].append(m)
                # Otherwise, create the appropriate R event
                else:
                    r_event = (NodeType.REARRANGEMENT, node, m)
                    event_graph[r_event] = [m]
                    event_graph[assignment].append(r_event)
    event_graph[node] = assignments

def create_o_events(node, G, event_graph):
    """
    Interpolate an Origin event node before each choice of location for an originating gene node
    """
    children = G[node]
    s_choice = children[1]
    root_syntenies = G[s_choice]
    event_graph[s_choice] = []
    for r in root_syntenies:
        origin_event = (NodeType.ORIGIN_EVENT, r)
        event_graph[s_choice].append(origin_event)
        event_graph[origin_event] = [r]
    event_graph[node] = children

def build_event_graph(G):
    """
    Build a graph that explicitly represents R and O events as nodes.
    This graph will be larger than the typical reconciliation graph,
    but is needed for computing medians w.r.t. the event distance.
    """
    event_graph = {}
    for node in graph_search_order(G):
        children = G[node]
        if node[0] is NodeType.LOCATION_MAPPING:
            if node[2] != "*":
                create_r_events(node, G, event_graph)
            # Copy the location assignment nodes
            # This makes it easier to case out copying the location assignments when
            # the syntenic location is not *, since the location assignment node does
            # not keep track of the syntenic location of the parent
            else:
                event_graph[node] = G[node]
                for child in G[node]:
                    event_graph[child] = G[child]
        elif node[0] is NodeType.ORIGIN:
            create_o_events(node, G, event_graph)
        # All location assignment and location list nodes are replaced in the previous cases
        # All other nodes will be faithfully kept
        elif not node[0] is NodeType.LOCATION_ASSIGNMENT and not node[0] is NodeType.LOCATION_LIST:
            event_graph[node] = children
    return event_graph
