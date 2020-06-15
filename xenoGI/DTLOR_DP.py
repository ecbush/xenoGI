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
# Edited by Ross Mawhorter in 2020 to reduce the running time per Rose Liu's formulation.

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

def delta(synteny1, synteny2, O, R):
    if synteny1==synteny2:
        return 0
    elif synteny1=="*":
        return O
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

def DP(hostTree, parasiteTree, phi, locus_map, D, T, L, Origin, R):
    """ Takes a hostTree, parasiteTree, tip mapping function phi, a locus_map, 
        and duplication cost (D), transfer cost (T), loss cost (L), 
        origin cost (O) and rearrange cost(R) and returns the DTLOR graph in the form of a dictionary, 
        as well as a the number of maximum parsimony reconciliations. The notation and 
        dynamic programming algorithm are explained in the tech report.
        Cospeciation is assumed to cost 0. """

    A = {}  # A, C, O, and bestSwitch are all defined in tech report
    C = {}
    O = {}
    eventsDict = {} # Dictionary to keep track of events that correspond to the min cost reconciliation 
    bestSwitch = {} 
    Minimums = {} # Dictionary to keep track of minimum reconciliation costs
    oBest = {} # Dictionary to keep track of the lowest costing events in O
    bestSwitchLocations = {} # Dictionary to keep track of switch locations
    allsynteny=list(locus_map.values())
    Allsynteny=copy.deepcopy(allsynteny)
    Allsynteny.append("*")
    #print("The dimensions is %d by %d by %d by %d"%(len(postorder(parasiteTree, "pTop")),len(Allsynteny), len(Allsynteny),len(postorder(hostTree, "hTop"))))
    for ep in postorder(parasiteTree, "pTop"):
        _,vp,ep1,ep2 = parasiteTree[ep]
        # is vp a tip?
        if ep1 == None: # then ep2 == None too and vp is a tip!
            vpIsATip = True
        else:
            vpIsATip = False
        for top_is_star in [True, False]:    #loop over all the possible locus 
            for l_bottom in valid_star(top_is_star, allsynteny):  #for start and end vertex of gene edge
                for eh in postorder(hostTree, "hTop"):
                    _,vh,eh1,eh2 = hostTree[eh]
                    eventsDict[(vp, vh, top_is_star, l_bottom)] = []
                    oBest[(vp, vh, top_is_star, l_bottom)] = []
                    # is vh a tip?
                    if eh1 == None: # then eh2 == None too and vh is a tip!
                        vhIsATip = True
                    else:
                        vhIsATip = False
                    # Compute A(ep, eh)
                    if vhIsATip:
                        if vpIsATip and phi[vp] == vh and locus_map[vp]==l_bottom:
                            A[(ep, eh, top_is_star, l_bottom)] = 0
                            # Contemporary event to be added to eventsDict
                            Amin = [["C", (None, None, None, None), (None, None, None, None)]] 
                        else: 
                            A[(ep, eh, top_is_star, l_bottom)] = Infinity
                            Amin = [Infinity]
                    else: #vh is not a tip
                        # Compute S and create event list to add to eventsDict
                        if not vpIsATip:
                            COepeh=Infinity
                            coMin = []
                            lowest_cost_spec=[]
                            for l1 in valid(l_bottom, allsynteny):
                                for l2 in valid(l_bottom, allsynteny):
                                    l_is_star = l_bottom == "*"
                                    synteny_cost = delta(l_bottom, l1, Origin, R) + delta(l_bottom, l2, Origin, R)
                                    # Add in the delta for handle branch 
                                    if ep=='pTop':
                                        synteny_cost += delta(top_is_star, l_bottom, Origin, R)

                                    #TODO these two lines take up >50% of time in the nested loop for locs
                                    co1=C[(ep1, eh1, l_is_star, l1)] + C[(ep2, eh2, l_is_star, l2)]
                                    co2=C[(ep1, eh2, l_is_star, l1)] + C[(ep2, eh1, l_is_star, l2)]
                                    
                                    CO_cost=min(co1,co2)
                                    CO_total=synteny_cost+CO_cost
                                    if CO_total<COepeh:
                                        COepeh=CO_total
                                        if CO_cost==co1:
                                            lowest_cost_spec=["S", (ep1, eh1, l_is_star, l1), \
                                            (ep2, eh2, l_is_star, l2)]
                                        else:
                                            lowest_cost_spec=["S", (ep1, eh2, l_is_star, l1), \
                                            (ep2, eh1, l_is_star, l2)]
                            if COepeh<Infinity:
                                
                                coMin.append(lowest_cost_spec)
                            else:
                                coMin = [Infinity]

                        else:
                            COepeh = Infinity
                            coMin = [Infinity]
                        # Compute L and create event list to add to eventsDict
                        lossMin = [] # List to keep track of lowest cost loss
                        loss_eh2=C[(ep, eh1, top_is_star, l_bottom)] #eh2 no longer has ep
                        loss_eh1=C[(ep, eh2, top_is_star, l_bottom)]
                        LOSSepeh = min(loss_eh2, loss_eh1)
                        #NOTE loss events record the mapping of parasite edge onto surviving host child
                        if LOSSepeh == loss_eh2: 
                            lossMin.append(["L", (vp, eh1, top_is_star, l_bottom), (None, None, None, None)])
                        if LOSSepeh == loss_eh1: 
                            lossMin.append(["L", (vp, eh2, top_is_star, l_bottom), (None, None, None, None)])
                        if ep!="pTop" and not top_is_star:
                            LOSSepeh += L 
                        # Determine which event occurs for A[(ep, eh)]
                        
                        A[(ep, eh, top_is_star, l_bottom)] = min(COepeh, LOSSepeh)     
                        # Record event occuring for A[(ep, eh)] as Amin
                        if COepeh < LOSSepeh:
                            Amin = coMin
                        elif LOSSepeh < COepeh: 
                            Amin = lossMin
                        else:
                            #lossMin, coMin are  list of lists
                            Amin = lossMin + coMin
                        

                    # Compute C(ep, eh,l_top, l_bottom)
                    #   First, compute D
                    if not vpIsATip:
                        DUPepeh=Infinity
                        # List to keep track of lowest cost duplication event
                        dupList=[Infinity]
                        for l1 in valid(l_bottom, allsynteny):
                            for l2 in valid(l_bottom, allsynteny):
                                l_is_star = l_bottom == "*"
                                #TODO: this line takes over half of time in this nested loop
                                dup_cost=delta(l_bottom, l1, Origin, R)+delta(l_bottom, l2, Origin, R) + \
                                        C[(ep1, eh, l_is_star, l1)]+C[(ep2, eh, l_is_star, l2)]
                                #only add the duplication cost if the bottom synteny on ep is not *
                                if l_bottom!="*":
                                    dup_cost+=D
                                if dup_cost<DUPepeh:
                                    DUPepeh=dup_cost
                                    dupList=["D", (ep1, vh, l_is_star, l1), (ep2, vh, l_is_star, l2)]
                    else:
                        DUPepeh = Infinity
                        dupList = [Infinity]
                   
                    #   Next, Compute T and create event list to add 
                    #   to eventsDict using bestSwitchLocations
                    if not vpIsATip and l_bottom!="*":
                        switchList = [] # List to keep track of lowest cost switch
                        SWITCHepeh=Infinity
                        #need to find all possible children syntenies
                        for l1 in valid(l_bottom, allsynteny):
                            for l2 in valid(l_bottom, allsynteny):
                                l_is_star = l_bottom == "*"
                                ep2_switch=C[(ep1, eh, l_is_star, l1)] + bestSwitch[(ep2, eh, l2)]
                                ep1_switch=C[(ep2, eh, l_is_star, l2)] + bestSwitch[(ep1, eh, l1)]
                                min_switch=min(ep2_switch, ep1_switch)
                                switch_cost = T+delta(l_bottom,l1, Origin, R)+delta(l_bottom,l2, Origin, R)+ min_switch
                                if switch_cost<SWITCHepeh:
                                    SWITCHepeh=switch_cost
                                    # if ep2 switching has the lowest cost
                                    if ep2_switch < ep1_switch:
                                        for location in bestSwitchLocations[(ep2,vh,l2)]:
                                            currentLoc = location[1] # Switch landing site
                                            switchList.append(["T", (ep1, vh, l1), (ep2, \
                                                currentLoc, l2)])
                                    # if ep1 switching has the lowest cost
                                    elif ep1_switch<ep2_switch: 
                                        for location in bestSwitchLocations[(ep1,vh, l1)]:
                                            currentLoc = location[1]
                                            switchList.append(["T", (ep2, vh, l2), \
                                                (ep1, currentLoc, l1)])
                                    # if ep1 switching has the same cost as ep2 switching
                                    else: 
                                        for location in bestSwitchLocations[(ep2, vh, l2)]:
                                            currentLoc = location[1]
                                            if currentLoc != None:
                                                switchList.append(["T", (ep1, vh, l1), (ep2, \
                                                currentLoc, l2)])
                                            else:
                                                switchList.append(["T", (ep1, vh, l1), \
                                                    (ep2, currentLoc, l2)])
                                   
                                        for location in bestSwitchLocations[(ep1,vh, l1)]:
                                            currentLoc = location[1]
                                            if currentLoc != None:
                                                switchList.append(["T", (ep2, vh, l2), \
                                                (ep1, currentLoc, l1)])
                                            else:
                                                switchList.append(["T", (ep1, vh, l1), \
                                                    (ep2, currentLoc, l2)])           
                        if switchList==[]:
                            switchList=[[Infinity]]

                    else:
                        SWITCHepeh = Infinity
                        switchList = [[Infinity]]
                    # Compute C[(ep, eh, l_top, l_bottom)] and add the event or events with that cost
                    # to the dictionary eventsDict
                    co_min=min(A[(ep, eh, top_is_star, l_bottom)], DUPepeh, SWITCHepeh)  
                    C[(ep, eh, top_is_star, l_bottom)] = co_min
                    Minimums[(vp, vh, top_is_star, l_bottom)] =co_min #min cost of reconciliation with this mapping and below
                    if co_min == DUPepeh:
                        eventsDict[(vp, vh, top_is_star, l_bottom)].append(dupList)  #duplist should just be a 1d list
                    if co_min == SWITCHepeh:
                        eventsDict[(vp, vh, top_is_star, l_bottom)].extend(switchList)     #switchList should be a list of lists
                    if co_min == A[(ep, eh, top_is_star, l_bottom)]:
                        eventsDict[(vp, vh, top_is_star, l_bottom)].extend(Amin)
                   
                    #do not allow top of gene tree handle to be an actual synteny
                    if Minimums[(vp, vh, top_is_star, l_bottom)] == Infinity or (ep=="pTop" and not top_is_star):
                        del Minimums[(vp, vh, top_is_star, l_bottom)]
                        del eventsDict[(vp, vh, top_is_star, l_bottom)]
                    # Compute O(ep, eh, l_top, l_bottom)
                    # Compute oBest[(vp, vh, l_top, l_bottom)], the source of O(ep, eh, l_top, l_bottom)
                    if vhIsATip: 
                        O[(ep, eh, top_is_star, l_bottom)] = C[(ep, eh, top_is_star, l_bottom)]  
                        oBest[(vp, vh, top_is_star, l_bottom)] = [(vp, vh, top_is_star, l_bottom)]              
                    else: 
                    
                        #finds Minimum Cost for O
                        O_list= [C[(ep, eh, top_is_star, l_bottom)], O[(ep, eh1, top_is_star, l_bottom)], O[(ep, eh2, top_is_star, l_bottom)]]
                        O_min=min(O_list)
                        O[(ep, eh, top_is_star, l_bottom)] = O_min     
                        #finds the minimum switch locations for O
                        oMin = O_list.index(O_min)
                        if oMin == 0:
                            oBest[(vp,vh, top_is_star, l_bottom)].append((vp, vh, top_is_star, l_bottom))
                        if oMin == 1:
                            oBest[(vp,vh, top_is_star, l_bottom)].extend(oBest[(vp, eh1, top_is_star, l_bottom)])
                        if oMin == 2:
                            oBest[(vp,vh, top_is_star, l_bottom)].extend(oBest[(vp, eh2, top_is_star, l_bottom)])
                # Compute bestSwitch values
                bestSwitch[(ep, "hTop", l_bottom)] = Infinity
                bestSwitchLocations[(vp, hostTree["hTop"][1], l_bottom)] = [(None,None, None, None)]
                for eh in preorder(hostTree, "hTop"):
                    _, vh, eh1, eh2 = hostTree[eh]

                    # is vh a tip?
                    if eh1 == None: # then eh2 == None too and vh is a tip!
                        vhIsATip = True
                    else:
                        vhIsATip = False
                    # find best place for a switch to occur (bestSwitch)
                    # and the location to which the edge switches (bestSwitchLocations)   
                    if eh1 != None and eh2 != None: # not a tip
                        bestSwitchLocations[(vp, eh1, l_bottom)] = []
                        bestSwitchLocations[(vp, eh2, l_bottom)] = []
                        ep_bestSwitch=bestSwitch[(ep, eh, l_bottom)]
                        O_eh2=O[(ep, eh2, top_is_star, l_bottom)]
                        O_eh1=O[(ep, eh1, top_is_star, l_bottom)]
                        bestSwitch[(ep, eh1, l_bottom)] = min(ep_bestSwitch,O_eh2)
                        bestSwitch[(ep, eh2, l_bottom)] = min(ep_bestSwitch,O_eh1)
                    
                        if bestSwitch[(ep, eh1, l_bottom)] == ep_bestSwitch and \
                        bestSwitchLocations[(vp, vh, l_bottom)] != [(None, None, None, None)]:
                            bestSwitchLocations[(vp, eh1, l_bottom)].extend\
                            (bestSwitchLocations[(vp, vh, l_bottom)])
                        if bestSwitch[(ep, eh1, l_bottom)] == O_eh2 and \
                        oBest[(vp, eh2, top_is_star, l_bottom)]!= [(None, None, None, None)]:
                            bestSwitchLocations[(vp, eh1, l_bottom)].extend\
                            (oBest[(vp, eh2, top_is_star, l_bottom)])
                        if bestSwitch[(ep, eh2, l_bottom)] == ep_bestSwitch and \
                        bestSwitchLocations[(vp, vh,  l_bottom)] != [(None, None, None, None)]:
                            bestSwitchLocations[(vp, eh2, l_bottom)].extend\
                            (bestSwitchLocations[(vp, vh, l_bottom)])
                        if bestSwitch[(ep, eh2, l_bottom)] == O_eh1 and \
                        oBest[(vp, eh1, top_is_star, l_bottom)]!=[(None, None, None, None)]:
                            bestSwitchLocations[(vp, eh2, l_bottom)].extend\
                            (oBest[(vp, eh1, top_is_star, l_bottom)])
    
    for key in bestSwitchLocations:
        if bestSwitchLocations[key][0] == (None, None, None, None):
            bestSwitchLocations[key] = bestSwitchLocations[key][1:]

    # Use findPath and findBestRoots to construct the DTLOR graph dictionary
    #this finds the optimal solutions (mappings)
    treeMin, min_cost = findBestRoots(parasiteTree, Minimums)
    #This picks a random MPR from the optimal ones
    MPR = findOneMPR(treeMin, eventsDict, {})  
    return MPR, min_cost


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
        if root != (None, None): #format (key, level)
            vertices = root[0] #format (vp, vh, l, l)
            if root[1] == 0:
                parentsDict[vertices] = ScoreDict[vertices]
            for n in range(len(DTLORDict[vertices])-1): #the last item is always the cost
                _,child1,child2,oldScore = DTLORDict[vertices][n]
                newDTLOR[vertices][n][3] = parentsDict[vertices] * \
                (1.0 * oldScore / ScoreDict[vertices])
                if child1!= (None, None, None, None):
                    if child1 in parentsDict:
                        parentsDict[child1] += newDTLOR[vertices][n][3]
                    else: 
                        parentsDict[child1] = newDTLOR[vertices][n][3] 
                if child2!=(None, None, None, None):
                    if child2 in parentsDict:
                        parentsDict[child2] += newDTLOR[vertices][n][3]
                    else: 
                        parentsDict[child2] = newDTLOR[vertices][n][3]
    normalize = newDTLOR[preOrder[-1][0]][0][-1]  #score of the top most event
    for key in newDTLOR:
        for event in newDTLOR[key][:-1]:
            event[-1] = event[-1]/normalize
    
    return newDTLOR, normalize

def findBestRoots(Parasite, MinimumDict):
    """Takes Parasite Tree and a dictionary of minimum reconciliation costs
    and returns a list of the minimum cost reconciliation tree roots and the associated minimum cost"""
    treeTops = []
    for key in MinimumDict:
        #pick out the set of keys in C (MinimumDict) that correspond to mapping 
        #the root of parasite to any place in host tree
        #and any locus combinations
        if key[0] == Parasite['pTop'][1]:
            treeTops.append(key)
    treeMin = []
    min_cost=min([MinimumDict[root] for root in treeTops])
    for pair in treeTops:
        #pick out the set of keys that gives the minimum of C (cost of optimal solution)
        #note there might be multiple equally optimal mappings
        if MinimumDict[pair] == min_cost:
            treeMin.append(pair)
    return treeMin,min_cost

def findPath(tupleList, eventDict, uniqueDict):
    """Takes as input tupleList, a list of minimum reconciliation cost mappings 
     of the form (vp, vh, l_top, l_bottom), eventDict, the dictionary of events and 
     children for each node, and uniqueDict, the dictionary of unique vertex mappings. 
     This returns the completed DTL graph as a Dictionary"""
     #findPath(treeMin, eventsDict, {})
    for mapping in tupleList:
        if not mapping in uniqueDict:
            uniqueDict[mapping] = eventDict[mapping]
        for event in eventDict[mapping][:-1]:
            for location in event:
                if type(location) is tuple and location != (None, None, None, None):
                    findPath([location], eventDict, uniqueDict)
    return uniqueDict

def findOneMPR(tupleList, eventDict, MPR):
    """
    Takes as input tupleList, a list of minimum reconciliation cost mappings 
     of the form (vp, vh, l_top, l_bottom), eventDict, the dictionary of events and 
     children for each node, and a dictionary to record one MPR, which is returned. 
    """
    mapping = choice(tupleList)
    event = choice(eventDict[mapping])  #randomly pick an event
    #print(event)
    MPR[mapping]=event  
    for location in event: #there should be eventType,and two child mapping nodes
        if type(location) is tuple and location != (None, None, None, None):
            findOneMPR([location], eventDict, MPR)
    return MPR

