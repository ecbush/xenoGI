# DTLOR_DP.py
# Nuo Liu, HMC 2019-2020 Thesis project
# Modified from Ran Libeskind-Hadas(June 2015)

# The basic DP algorithm for reconciling pairs of trees
# Altered and expanded by Carter Slocum and Annalise Schweickart
# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop".

# Edited by Annalise Schweickart and Carter Slocum, July 2015 to return
# the DTL reconciliation graph that uses frequency scoring, as well as the
# number of reconciliations of the host and parasite trees

# Modified by Ross Mawhorter in 2020 to reduce the running time per Rose Liu's formulation
# and apply the new model of how DTLOR should work.

# import newickFormatReader
from .Greedy import *
import copy
import sys, glob, os
import time
from random import choice

def valid_star(is_top_star, allsynteny):
    """
    returns the set of valid syntenies for descendants
    """
    if is_top_star:
        new_locs = copy.deepcopy(allsynteny)
        new_locs.append("*")
        return new_locs
    else:
        return allsynteny

def valid(synteny, allsynteny):
    """
    returns the set of valid syntenies for the descendents
    """
    is_top_star = (synteny == "*")
    return valid_star(is_top_star, allsynteny)

def delta_cost(synteny1, synteny2, O, R):
    if synteny1==synteny2:
        return 0
    elif synteny1=="*":
        return O
    else:
        return R

# New delta function that does not use O since
# the O costs are calculated separately
def delta_r(synteny1, synteny2, R):
    if synteny1 == synteny2:
        return 0
    else:
        return R

Infinity = float('inf')

def preorder(tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns a 
    list of the edges in that tree in preorder (high edges to low edges)"""

    value = tree[rootEdgeName]
    _,_,leftChildEdgeName,rightChildEdgeName = value

    # base case
    if leftChildEdgeName == None: # then rightChildEdgeName == None also
        return [rootEdgeName]
    else:
        return [rootEdgeName] + \
                preorder(tree, leftChildEdgeName) + \
                preorder(tree, rightChildEdgeName)

def postorder(tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns a 
    list of the edges in that tree in postorder (low edges to high edges)"""

    value = tree[rootEdgeName]
    _,_,leftChildEdgeName,rightChildEdgeName = value
    # base case
    if leftChildEdgeName == None: # then rightChildEdgeName == None also
        return [rootEdgeName]
    else:
        return postorder(tree, leftChildEdgeName) + \
               postorder(tree, rightChildEdgeName) + \
               [rootEdgeName]

def check_tip(vh, eh1, eh2):
    if eh1 is None and eh2 is None:
        return True
    elif eh1 is not None and eh2 is not None:
        return False
    else:
        assert False, "Species node with one child: {}".format(vh)

#TODO: rename
def find_min_events(events_list):
    """
    Combines the events list into a single events list
    Each element of events_list is a tuple of (cost, [events])
    Which holds the DP entry for some part of the DP and the associated events
    This combines them into the entry of lowest cost, combining
    all events that have lowest cost
    """
    cost = Infinity
    events = []
    for c,e in events_list:
        if cost > c:
            cost = c
            events = e
        elif cost == c:
            events.extend(e)
    return (cost, events)


#TODO: use this instead of list comprehensions
# To do this, replace the list comprehensions in the current code with this
# then use the old find_min_events to combine the output once all the outputs
# are in the (cost, [event]) form
def find_min_events_alt(elements, cost_computer, event_computer):
    """
    Faster alternative to find_min_events.
    Does not compute events for a given cost unless that cost is optimal.
    Also should use way less memory than find_min_events.
    For example, if the elements are [(species_node, location)]
    Then cost_computer should be type:
      (species_node, location) -> float
    and event_computer should be type:
      (species_node, location) -> event
    """
    min_cost = Infinity
    min_events = []
    costs = map(lambda e: cost_computer(*e), elements)
    for element, cost in zip(elements, costs):
        if cost < min_cost:
            min_cost = c
            min_events = event_computer(*element)
        elif cost == min_cost:
            min_events.extend(event_computer(*element))
    return (min_cost, min_events)

#TODO: optimization where valid returns only the synteny locations for the clade UNDER a given gene node
#TODO: Refactor other code to reflect the tree representation change:
# ("C", None, None) versus ("C", (None, None, None, None), (None, None, None, None))
# and ("L", (...), None) versus ("L", (...), (None, None, None, None)

def DP(hostTree, parasiteTree, phi, locus_map, D, T, L, Origin, R):
    """ Takes a hostTree, parasiteTree, tip mapping function phi, a locus_map, 
        and duplication cost (D), transfer cost (T), loss cost (L), 
        origin cost (O) and rearrange cost(R) and returns the an MPR as a dict.
        The notation and dynamic programming algorithm are explained in the tech report.
        Cospeciation is assumed to cost 0. """

    A = {}  # A, C, O, and best_switch are all defined in tech report
    C = {}
    O = {}
    eventsDict = {} # Dictionary to keep track of events that correspond to the min cost reconciliation 
    best_switch = {} 
    allsynteny=list(locus_map.values())
    # Capture R for ease of use - it never changes
    delta = lambda s1, s2: delta_r(s1, s2, R)
    #print("The dimensions is %d by %d by %d by %d"%(len(postorder(parasiteTree, "pTop")),len(Allsynteny), len(Allsynteny),len(postorder(hostTree, "hTop"))))
    for ep in postorder(parasiteTree, "pTop"):
        _,vp,ep1,ep2 = parasiteTree[ep]
        vp_is_a_tip = check_tip(vp, ep1, ep2)
        for lp in allsynteny:  # The location of ep at the bottom of the branch above ep
            for eh in postorder(hostTree, "hTop"):
                _,vh,eh1,eh2 = hostTree[eh]
                eventsDict[(vp, vh, lp)] = []
                vh_is_a_tip = check_tip(vh, eh1, eh2)
                # Compute A[(ep, eh, lp)]
                if vh_is_a_tip:
                    if vp_is_a_tip and phi[vp] == vh and locus_map[vp]==lp:
                        A[(ep, eh, lp)] = (0, [("C", None, None)])
                    else: 
                        A[(ep, eh, lp)] = (Infinity, [])
                else: # vh is not a tip
                    # Compute cospeciation events
                    if not vp_is_a_tip:
                        def get_cospeciations(l1, l2):
                            synteny_cost = delta(lp, l1) + delta(lp, l2)
                            co1 = (synteny_cost + C[(ep1, eh1, l1)][0] + C[(ep2, eh2, l2)][0], \
                                    [("S", (ep1, eh1, l1), (ep2, eh2, l2))])
                            co2 = (synteny_cost + C[(ep1, eh2, l1)][0] + C[(ep2, eh1, l2)][0], \
                                    [("S", (ep1, eh2, l1), (ep2, eh1, l2))])
                            return find_min_events([co1, co2])
                        cospeciation_list = [get_cospeciations(l1, l2) for l1 in allsynteny for l2 in allsynteny]
                        cospeciations = find_min_events(cospeciation_list)
                    else:
                        cospeciations = (Infinity, [])
                    # Compute loss events
                    # eh1 is the branch where ep is lost
                    loss_eh1 = (C[(ep, eh2, lp)][0] + L, [("L", (vp, eh2, lp), None)])
                    # eh2 is the branch where ep is lost
                    loss_eh2 = (C[(ep, eh1, lp)][0] + L, [("L", (vp, eh1, lp), None)])
                    losses = find_min_events([loss_eh1, loss_eh2])

                    # Determine which event occurs for A[(ep, eh, lp)]
                    A[(ep, eh, lp)] = find_min_events([cospeciations, losses])

                # Compute C[(ep, eh,l_top, lp)]
                # First, compute duplications
                if not vp_is_a_tip:
                    DUPepeh=Infinity
                    # List to keep track of lowest cost duplication event
                    dupList=[Infinity]
                    def get_duplication(l1, l2):
                        synteny_cost = delta(lp, l1) + delta(lp, l2)
                        dup_cost = C[(ep1, eh, l1)][0] + C[(ep2, eh, l2)][0] + D
                        dup_event = ("D", (ep1, vh, l1), (ep2, vh, l2))
                        return (dup_cost, [dup_event])
                    dup_list = [get_duplication(l1, l2) for l1 in allsynteny for l2 in allsynteny]
                    duplications = find_min_events(dup_list)
                else:
                    duplications = (Infinity, [])
               
                # Compute transfer table
                if not vp_is_a_tip:
                    #TODO: the transferred child keeps the same synteny?
                    # Over all possible children syntenies
                    def get_transfer(new_l):
                        synteny_cost = delta(lp, new_l)
                        # Cost to transfer ep2
                        # Transferred child (ep2) keeps the same syntenic location (lp)
                        ep2_cost, ep2_locations = best_switch[(ep2, eh, lp)]
                        ep2_switch_cost = T + synteny_cost + C[(ep1, eh, new_l)][0] + ep2_cost
                        ep2_switch_events = [("T", (ep1, vh, new_l), (ep2, location[1], lp)) for location in \
                                ep2_locations]
                        ep2_switch = (ep2_switch_cost, ep2_switch_events)
                        # Cost to transfer ep1
                        # Now ep1 is being transferred and keeps lp
                        ep1_cost, ep1_locations = best_switch[(ep1, eh, lp)]
                        ep1_switch_cost = T + synteny_cost + C[(ep2, eh, new_l)][0] + ep1_cost
                        ep1_switch_events = [("T", (ep2, vh, new_l), (ep1, location[1], lp)) for location in \
                                ep1_locations]
                        ep1_switch = (ep1_switch_cost, ep1_switch_events)
                        return find_min_events([ep2_switch, ep1_switch])
                    transfer_list = [get_transfer(l1) for l1 in allsynteny]
                    transfers = find_min_events(transfer_list)
                else:
                    transfers = (Infinity, [])

                # Compute C[(ep, eh, l_top, lp)] and add the event or events with that cost
                # to the dictionary eventsDict
                C[(ep, eh, lp)] = \
                        find_min_events([A[(ep, eh, lp)], duplications, transfers])
                # The root must factor in the cost of getting a syntenic location
                if ep == "pTop":
                    old = C[(ep, eh, lp)]
                    C[(ep, eh, lp)] = (old[0] + Origin, old[1])

                # Compute O[(ep, eh, lp)]
                # O is mapping_node -> (cost, [mapping_node])
                if vh_is_a_tip: 
                    O[(ep, eh, lp)] = (C[(ep, eh, lp)][0], [(vp, vh, lp)])
                else: 
                    O_c = (C[(ep, eh, lp)][0], [(vp, vh, lp)])
                    O_eh1 = (O[(ep, eh1, lp)])
                    O_eh2 = (O[(ep, eh2, lp)])
                    O[(ep, eh, lp)] = find_min_events([O_c, O_eh1, O_eh2])

            # Compute best_switch values for the children
            best_switch[(ep, "hTop", lp)] = (Infinity, [])
            for eh in preorder(hostTree, "hTop"):
                _, vh, eh1, eh2 = hostTree[eh]
                # Find the best switches and switch locations for the children of vh
                # Don't set best_switch for nonexistent children
                if not check_tip(vh, eh1, eh2):
                    ep_best_switch = best_switch[(ep, eh, lp)]
                    O_eh2=O[(ep, eh2, lp)]
                    O_eh1=O[(ep, eh1, lp)]
                    best_switch[(ep, eh1, lp)] = find_min_events([ep_best_switch, O_eh2])
                    best_switch[(ep, eh2, lp)] = find_min_events([ep_best_switch, O_eh1])
        # Compute the cost of not giving a syntenic location
        # Tip must have a syntenic location
        if vp_is_a_tip:
            C[(vp, "hTop", "*")] = (Infinity, [])
        else:
            def get_single_null(eh, l):
                # Left child stays null
                l_map = (ep1, "hTop", "*")
                r_map = (ep2, eh, l)
                left_null_cost = C[l_map][0] + C[r_map][0] + Origin
                left_null_event = ("N", l_map, r_map)
                left_null = (left_null_cost, [left_null_event])
                # Right child stays null
                l_map = (ep1, eh, l)
                r_map = (ep2, "hTop", "*")
                right_null_cost = C[l_map][0] + C[r_map][0] + Origin
                right_null_event = ("N", l_map, r_map)
                right_null = (right_null_cost, [right_null_event])
                return find_min_events([left_null, right_null])
            single_null_list = [get_single_null(eh, l) \
                    for eh in postorder(hostTree, "hTop") for l in allsynteny]
            single_null = find_min_events(single_null_list)

            # Neither child gets a synteny
            l_map = (ep1, "hTop", "*")
            r_map = (ep2, "hTop", "*")
            both_null_cost = C[l_map][0] + C[r_map][0]
            both_null_event = ("N", l_map, r_map)
            both_null = (both_null_cost, [both_null_event])

            def get_neither_child_null(eh1, l1, eh2, l2):
                l_map = (ep1, eh1, l1)
                r_map = (ep2, eh2, l2)
                neither_null_cost = C[l_map][0] + C[r_map][0] + 2 * Origin
                neither_null_event = ("N", l_map, r_map)
                return (neither_null_cost, [neither_null_event])

            neither_null_list = [get_neither_child_null(eh1, l1, eh2, l2) \
                    for eh1 in postorder(hostTree, "hTop") for eh2 in postorder(hostTree, "hTop") \
                    for l1 in allsynteny for l2 in allsynteny]
            neither_null = find_min_events(neither_null_list)

            C[(ep, "hTop", "*")] = find_min_events([single_null, both_null, neither_null])

    # Cost for assigning the root a syntenic location
    root_not_null_list = [C[("pTop", eh, l)] for eh in postorder(hostTree, "hTop") for l in allsynteny]
    # Cost for not assigning a syntenic location
    root_null = C[("pTop", "hTop", "*")]
    root_list = root_not_null_list + [root_null]
    min_cost, _ = find_min_events(root_list)

    # Find the mapping nodes involving pTop of minimum cost
    best_roots = [m for m,c in C.items() if m[0] == "pTop" and c[0] == min_cost]

    # This picks a random MPR from the optimal ones
    MPR = find_MPR(best_roots, C)
    return MPR, min_cost

def find_MPR(best_roots, C):
    """
    Find a single MPR for the given C dict and best_roots.
    """
    return find_MPR_helper(best_roots, C, {})

def find_MPR_helper(nodes, C, MPR):
    """
    Recursively find a single MPR using C. Does the work for find_MPR.
    """
    mapping = choice(nodes)
    event = choice(C[mapping][1])
    MPR[mapping] = event
    e_type, e_left, e_right = event
    if e_left is not None:
        _ = find_MPR_helper([e_left], C, MPR)
    if e_right is not None:
        _ = find_MPR_helper([e_right], C, MPR)
    return MPR

def MPR_graph(best_roots, C):
    """
    Find the MPR graph for a given C dict and best_roots.
    """
    return MPR_graph_helper(best_roots, C, {})

def MPR_graph_helper(nodes, C, G):
    """
    Recursively create the entire MPR graph. Does the work for MPR_graph.
    """
    for mapping in nodes:
        events = c[mapping][1]
        G[mapping] = events
        for e_type, e_left, e_right in events:
            MPR_graph_helper([e_left], C, G)
            MPR_graph_helper([e_right], C, G)
    return G

def preorderDTLORsort(DTLOR, ParasiteRoot):
    """This takes in a DTL reconciliation graph and parasite root and returns 
    a sorted list, orderedKeysL, that is ordered by level from largest to 
    smallest, where level 0 is the root and the highest level has tips."""

    keysL = orderDTLOR(DTLOR, ParasiteRoot)
    uniqueKeysL = sortHelper(DTLOR, keysL)
    orderedKeysL = []
    levelCounter = 0
    while len(orderedKeysL) < len(keysL):
        for mapping in keysL:
            if mapping[-1] == levelCounter:
                orderedKeysL = orderedKeysL + [mapping]
        levelCounter += 1
    
    # lastLevel = orderedKeysL[-1][1]
    return orderedKeysL

def addScores(treeMin, DTLORDict, ScoreDict):
    """Takes the list of reconciliation roots, the DTLOR reconciliation graph, 
    a dictionary of parent nodes, and a dictionary of score values, and 
    returns the DTLOR with the normalized frequency scores calculated."""
    newDTLOR = copy.deepcopy(DTLORDict)
    parentsDict = {}
    ParasiteRoot=treeMin[0][0]
    preOrder = preorderDTLORsort(DTLORDict, ParasiteRoot)
    for root in preOrder:
        #TODO: when would this be true?
        if root != (None, None): #format (key, level)
            vertices = root[0] #format (vp, vh, l, l)
            if root[1] == 0:
                parentsDict[vertices] = ScoreDict[vertices]
            for n in range(len(DTLORDict[vertices])-1): #the last item is always the cost
                _,child1,child2,oldScore = DTLORDict[vertices][n]
                newDTLOR[vertices][n][3] = parentsDict[vertices] * \
                (1.0 * oldScore / ScoreDict[vertices])
                if child1 is not None:
                    if child1 in parentsDict:
                        parentsDict[child1] += newDTLOR[vertices][n][3]
                    else: 
                        parentsDict[child1] = newDTLOR[vertices][n][3] 
                if child2 is not None:
                    if child2 in parentsDict:
                        parentsDict[child2] += newDTLOR[vertices][n][3]
                    else: 
                        parentsDict[child2] = newDTLOR[vertices][n][3]
    normalize = newDTLOR[preOrder[-1][0]][0][-1]  #score of the top most event
    for key in newDTLOR:
        for event in newDTLOR[key][:-1]:
            event[-1] = event[-1]/normalize
    
    return newDTLOR, normalize

