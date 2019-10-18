# DTLOR_DP.py
# Nuo Liu
# Modified from Ran Libeskind-Hadas, June 2015

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

# import newickFormatReader
import Greedy
import copy
from treeParser import parseTreeForDP
from  .familiesDTLOR import *
from .Family import *


def allSynteny(familiesO):
    #synteny numbering start from 0
    return list(range(familiesO.getNumLocusFamilies()))

def valid(synteny, allsynteny):
    """
    returns the set of valid syntenies for the descendents
    """
    if synteny=="*": 
        return allsynteny.append("*")
    else: 
        return allsynteny

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

def DP(hostTree, parasiteTree, phi, locus_map, allsynteny, D, T, L, O, R):
    """ Takes a hostTree, parasiteTree, tip mapping function phi, a locus_map, all unique
        syntenic locations, and duplication cost (D), transfer cost (T), loss cost (L), 
        origin cost (O) and rearrange cost(R) and returns the DTLOR graph in the form of a dictionary, 
        as well as a the number of maximum parsimony reconciliations. The notation and 
        dynamic programming algorithm are explained in the tech report.
        Cospeciation is assumed to cost 0. """
    #preprocess the trees into the correct format
    hostTree=parseTreeForDP(hostTree,parasite=False)
    parasiteTree=parseTreeForDP(parasiteTree,parasite=True)


    A = {}  # A, C, O, and bestSwitch are all defined in tech report
    C = {}
    O = {}
    eventsDict = {} # Dictionary to keep track of events, children, and scores
    bestSwitch = {} 
    Minimums = {} # Dictionary to keep track of minimum reconciliation costs
    oBest = {} # Dictionary to keep track of the lowest costing events in O
    bestSwitchLocations = {} # Dictionary to keep track of switch locations
    Score = {} # Dictionary to calculate the frequency scoring of each event
    all_loc=allsynteny.append("*")
    # Following logic taken from tech report
    for ep in postorder(parasiteTree, "pTop"):
        for eh in postorder(hostTree, "hTop"):
            for l_top in all_loc:    #loop over all the possible locus 
                for l_bottom in all_loc:  #for start and end vertex of gene edge
            
                    _,vp,ep1,ep2 = parasiteTree[ep]
                    _,vh,eh1,eh2 = hostTree[eh]
                    eventsDict[(vp, vh, l_top, l_bottom)] = []
                    oBest[(vp, vh, l_top, l_bottom)] = []
                    # is vp a tip?
                    if ep1 == None: # then ep2 == None too and vp is a tip!
                        vpIsATip = True
                        pChild1 = None
                        pChild2 = None
                    else:
                        vpIsATip = False
                        #TODO figure out why there is another indexing
                        pChild1 = parasiteTree[ep][2][1]
                        pChild2 = parasiteTree[ep][3][1]

                    # is vh a tip?
                    if eh1 == None: # then eh2 == None too and vh is a tip!
                        vhIsATip = True
                        hChild1 = None
                        hChild2 = None
                    else:
                        vhIsATip = False
                        hChild1 = hostTree[eh][2][1]
                        hChild2 = hostTree[eh][3][1]
                        
                    # Compute A(ep, eh)

                    if vhIsATip:
                        if vpIsATip and phi[vp] == vh and locus_map[vp]==l_bottom:
                            A[(ep, eh, l_top, l_bottom)] = 0
                            # Contemporary event to be added to eventsDict
                            #TODO: figure out what this should change to ????
                            Amin = [["C", (None, None), (None, None), 1]] 
                            Score[(vp, vh, l_top, l_bottom)] = 1.0
                        else: 
                            Score[(vp, vh,l _top, l_bottom)] = Infinity
                            A[(ep, eh, l_top, l_bottom)] = Infinity
                            Amin = [Infinity]
                    else: #vh is not a tip
                        # Compute S and create event list to add to eventsDict
                        if not vpIsATip:
                            COepeh=Infinity
                            coMin = []
                            lowest_cost_spec=[]
                            for l1 in valid(l_bottom):
                                for l2 in valid(l_bottom):
                                    synteny_cost=delta(l_bottom,l1)+delta(l_bottom,l2)
                                    CO_cost=min(C[(ep1, eh1,l_bottom, l1)] + C[(ep2, eh2, l_bottom, l2)], \
                                        C[(ep1, eh2, l_bottom, l1)] + C[(ep2, eh1, l_bottom, l2)])
                                    CO_total=synteny_cost+CO_cost
                                    if CO_total<COepeh:
                                        COepeh=CO_total
                                        if CO_cost==C[(ep1, eh1,l_bottom, l1)] + C[(ep2, eh2, l_bottom, l2)]:
                                            lowest_cost_spec=["S", (pChild2, hChild1,l_bottom, l1), \
                                            (pChild1, hChild2, l_bottom, l2), (Score[(pChild2, hChild1, l_bottom, l1)] * \
                                            Score[(pChild1, hChild2, l_bottom, l2)])]
                                        if COepeh == C[(ep1, eh2, l_bottom, l1)] + C[(ep2, eh1, l_bottom, l2)]:
                                            lowest_cost_spec=["S", (pChild1, hChild1,l_bottom,l1), \
                                            (pChild2, hChild2, l_bottom, l2),(Score[(pChild1, hChild1, l_bottom, l1)]\
                                             * Score[(pChild2, hChild2, l_bottom, l2)])]
                            coMin.append(lowest_cost_spec)
                        else:
                            COepeh = Infinity
                            coMin = [Infinity]
                            Score[(vp, vh, l_top, l_bottom)] = Infinity
                        # Compute L and create event list to add to eventsDict
                        lossMin = [] # List to keep track of lowest cost loss
                        LOSSepeh = min(C[(ep, eh1, l_top, l_bottom)], C[(ep, eh2, l_top, l_bottom)])
                            #NOTE loss events record the mapping of parasite edge onto surviving host child
                            if LOSSepeh == C[(ep, eh1, l_top, l_bottom)]: lossMin.append(\
                                ["L", (vp, hChild1, l_top, l_bottom), (None, None, None, None), \
                                Score[(vp, hChild1, l_top, l_bottom)]])
                            if LOSSepeh == C[(ep, eh2, l_top, l_bottom)]: lossMin.append(\
                                ["L", (vp, hChild2,l_top, l_bottom), (None, None, None, None), Score[(vp, hChild2,l_top, l_bottom)]])
                        if ep!="pTop":
                            LOSSepeh = L + min(C[(ep, eh1, l_top, l_bottom)], C[(ep, eh2, l_top, l_bottom)])

                        # Determine which event occurs for A[(ep, eh)]
                        A[(ep, eh, l_top, l_bottom)] = min(COepeh, LOSSepeh)
                        # Record event occuring for A[(ep, eh)] as Amin
                        if COepeh < LOSSepeh:
                            Amin = coMin
                        elif LOSSepeh < COepeh: 
                            Amin = lossMin
                        else:
                            Amin = lossMin + coMin

                    # Compute C(ep, eh,l_top, l_bottom)
                    #   First, compute D
                    if not vpIsATip:
                        DUPepeh=Infinity
                        # List to keep track of lowest cost duplication event
                        dupList=[]
                        for l1 in valid(l_bottom):
                            for l2 in valid(l_bottom):
                                dup_cost=D+delta(l_bottom,l1)+delta(l_bottom,l2)+ \
                                    C[(ep1, eh, l_bottom, l1)]+C[(ep2, eh, l_bottom, l2)]
                               
                                if dup_cost<DUPepeh:
                                    DUPepeh=dup_cost
                                    dupList=["D", (pChild1, vh, l_bottom, l1), (pChild2, vh, l_bottom, l2), \
                                            (Score[(pChild1, vh, l_bottom, l1)] * Score[(pChild2, vh,l_bottom, l2)])]
                    else:
                        DUPepeh = Infinity
                        dupList = [Infinity]
                   
                    #   Next, Compute T and create event list to add 
                    #   to eventsDict using bestSwitchLocations
                    if not vpIsATip and l_bottom!="*":
                        switchList = [] # List to keep track of lowest cost switch
                        SWITCHepeh=Infinity
                        for l1 in valid(l_bottom):
                            for l2 in valid(l_bottom):
                                switch_cost = T + delta(l_bottom,l1)+delta(l_bottom,l2)+ \
                                    min(C[(ep1, eh, l_bottom, l1)] + bestSwitch[(ep2, eh, l_bottom, l2)],C[(ep2, eh, l_bottom, l2)] + bestSwitch[(ep1, eh, l_bottom, l1)]) 
                                if switch_cost<SWITCHepeh:
                                    SWITCHepeh=switch_cost
                                    # if ep2 switching has the lowest cost
                                    if (C[(ep1, eh, l_bottom, l1)] + bestSwitch[(ep2, eh, l_bottom, l2)]) < \
                                        (C[(ep2, eh, l_bottom, l2)] + bestSwitch[(ep1, eh, l_bottom, l1)]):
                                        for location in bestSwitchLocations[(pChild2,vh,l_bottom,l2)]:
                                            currentLoc = location[1] # Switch landing site
                                            if currentLoc == None: # Switches to a leaf
                                                Score[(pChild1, currentLoc, l_bottom, l1)] = Infinity
                                                Score[(pChild2, currentLoc, l_bottom, l2)] = Infinity
                                            switchList.append(["T", (pChild1, vh, l_bottom, l1), (pChild2, \
                                                currentLoc, l_bottom, l2), (Score[(pChild1, vh, l_bottom, l1)] * \
                                                Score[(pChild2, currentLoc, l_bottom, l2)])])
                                    # if ep1 switching has the lowest cost
                                    elif (C[(ep2, eh, l_bottom, l2)] + bestSwitch[(ep1, eh, l_bottom, l1)]) < \
                                        (C[(ep1, eh, l_bottom, l1)] + bestSwitch[(ep2, eh, l_bottom, l2)]): 
                                        for location in bestSwitchLocations[(pChild1,vh, l_bottom,l1)]:
                                            currentLoc = location[1]
                                            if currentLoc == None:
                                                Score[(pChild1, currentLoc, l_bottom, l1)] = Infinity
                                                Score[(pChild2, currentLoc, l_bottom, l2)] = Infinity
                                            switchList.append(["T", (pChild2, vh, l_bottom, l2), \
                                                (pChild1, currentLoc, l_bottom, l1), (Score[(pChild2, vh, l_bottom, l2)] * \
                                                    Score[(pChild1, currentLoc, l_bottom, l1)])])
                                    # if ep1 switching has the same cost as ep2 switching
                                    else: 
                                        for location in bestSwitchLocations[(pChild2,vh,l_bottom,l2)]:
                                            currentLoc = location[1]
                                            if currentLoc != None:
                                                switchList.append(["T", (pChild1, vh, l_bottom, l1), (pChild2, \
                                                currentLoc, l_bottom, l2), (Score[(pChild1, vh, l_bottom, l1)] * \
                                                Score[(pChild2, currentLoc, l_bottom, l2)])])
                                            else:
                                                switchList.append(["T", (pChild1, vh, l_bottom, l1), \
                                                    (pChild2, currentLoc, l_bottom, l2), Infinity])
                                        for location in bestSwitchLocations[(pChild1,vh, l_bottom,l1)]:
                                            currentLoc = location[1]
                                            if currentLoc != None:
                                                switchList.append(["T", (pChild2, vh, l_bottom, l2), \
                                                (pChild1, currentLoc, l_bottom, l1), (Score[(pChild2, vh, l_bottom, l2)] * \
                                                    Score[(pChild1, currentLoc, l_bottom, l1)])])
                                            else:
                                                switchList.append(["T", (pChild1, vh, l_bottom, l1), \
                                                    (pChild2, currentLoc, l_bottom, l2), Infinity])

                    else:
                        SWITCHepeh = Infinity
                        switchList = [Infinity]
                    # Compute C[(ep, eh, l_top, l_bottom)] and add the event or events with that cost
                    # to the dictionary eventsDict
                    C[(ep, eh, l_top, l_bottom)] = min(A[(ep, eh, l_top, l_bottom)], DUPepeh, SWITCHepeh)
                    Minimums[(vp, vh, l_top, l_bottom)] = C[(ep, eh, l_top, l_bottom)]
                    if C[(ep, eh, l_top, l_bottom)] == DUPepeh:
                        eventsDict[(vp, vh, l_top, l_bottom)].append(dupList)
                    if C[(ep, eh, l_top, l_bottom)] == SWITCHepeh:
                        eventsDict[(vp, vh, l_top, l_bottom)].extend(switchList)
                    if C[(ep, eh, l_top, l_bottom)] == A[(ep, eh, l_top, l_bottom)]:
                        eventsDict[(vp, vh, l_top, l_bottom)].extend(Amin)
                    for key in eventsDict:
                        mapScore = 0 # initialize frequency scoring for each event
                        for event in eventsDict[key]:
                            if type(event) is list:
                                mapScore += event[-1]
                        Score[key] = mapScore
                    if Minimums[(vp, vh, l_top, l_bottom)] == Infinity:
                        del Minimums[(vp, vh, l_top, l_bottom)]
                        del eventsDict[(vp, vh, l_top, l_bottom)]
                    # Compute O(ep, eh, l_top, l_bottom)
                    # Compute oBest[(vp, vh, l_top, l_bottom)], the source of O(ep, eh, l_top, l_bottom)
                    if vhIsATip: 
                        O[(ep, eh, l_top, l_bottom)] = C[(ep, eh, l_top, l_bottom)]  
                        oBest[(vp, vh, l_top, l_bottom)] = [(vp, vh, l_top, l_bottom)]              
                    else: 
                    
                        #finds Minimum Cost for O
                        O[(ep, eh, l_top, l_bottom)] = min(C[(ep, eh, l_top, l_bottom)], O[(ep, eh1, l_top, l_bottom)], O[(ep, eh2, l_top, l_bottom)])
                        #finds the minimum switch locations for O
                        oMin = [C[(ep, eh, l_top, l_bottom)], O[(ep, eh1, l_top, l_bottom)], O[(ep, eh2, l_top, l_bottom)]].index\
                                                                                                (O[(ep, eh, l_top, l_bottom)])
                        if oMin == 0:
                            oBest[(vp,vh, l_top, l_bottom)].append((vp, vh, l_top, l_bottom))
                        if oMin == 1:
                            oBest[(vp,vh, l_top, l_bottom)].extend(oBest[(vp, hChild1, l_top, l_bottom)])
                        if oMin == 2:
                            oBest[(vp,vh, l_top, l_bottom)].extend(oBest[(vp, hChild2, l_top, l_bottom)])

                    
                        
                    # Compute bestSwitch values
                    bestSwitch[(ep, "hTop", l_top, l_bottom)] = Infinity
                    bestSwitchLocations[(vp, hostTree["hTop"][1], l_top, l_bottom)] = [(None,None, None, None)]
                    for eh in preorder(hostTree, "hTop"):
                        _, vp, ep1, ep2 = parasiteTree[ep]
                        _, vh, eh1, eh2 = hostTree[eh]

                        #is vp a tip?
                        if ep1 == None:
                            vpIsATip = True
                            pChild1 = None
                            pChild2 = None
                        else:
                            vpIsATip = False
                            pChild1 = parasiteTree[ep][2][1]
                            pChild2 = parasiteTree[ep][3][1]

                        # is vh a tip?
                        if eh1 == None: # then eh2 == None too and vh is a tip!
                            vhIsATip = True
                            hChild1 = None
                            hChild2 = None
                        else:
                            vhIsATip = False
                            hChild1 = hostTree[eh][2][1]
                            hChild2 = hostTree[eh][3][1]
                        # find best place for a switch to occur (bestSwitch)
                        # and the location to which the edge switches (bestSwitchLocations)   
                        if eh1 != None and eh2 != None: # not a tip
                            bestSwitchLocations[(vp, hChild1, l_top, l_bottom)] = []
                            bestSwitchLocations[(vp, hChild2, l_top, l_bottom)] = []
                            bestSwitch[(ep, eh1, l_top, l_bottom)] = min(bestSwitch[(ep, eh, l_top, l_bottom)],\
                                                                        O[(ep, eh2, l_top, l_bottom)])
                            bestSwitch[(ep, eh2, l_top, l_bottom)] = min(bestSwitch[(ep, eh, l_top, l_bottom)],\
                                                                        O[(ep, eh1, l_top, l_bottom)])
                            if bestSwitch[(ep, eh1, l_top, l_bottom)] == bestSwitch[(ep, eh, l_top, l_bottom)] and \
                            bestSwitchLocations[(vp, vh, l_top, l_bottom)] != [(None, None, None, None)]:
                                bestSwitchLocations[(vp, hChild1, l_top, l_bottom)].extend\
                                (bestSwitchLocations[(vp, vh, l_top, l_bottom)])
                            if bestSwitch[(ep, eh1, l_top, l_bottom)] == O[(ep, eh2, l_top, l_bottom)] and \
                            oBest[(vp, hChild2, l_top, l_bottom)]!= [(None, None, None, None)]:
                                bestSwitchLocations[(vp, hChild1, l_top, l_bottom)].extend\
                                (oBest[(vp, hChild2, l_top, l_bottom)])
                            if bestSwitch[(ep, eh2, l_top, l_bottom)] == bestSwitch[(ep, eh, l_top, l_bottom)] and \
                            bestSwitchLocations[(vp, vh, l_top, l_bottom)] != [(None, None, None, None)]:
                                bestSwitchLocations[(vp, hChild2, l_top, l_bottom)].extend\
                                (bestSwitchLocations[(vp, vh, l_top, l_bottom)])
                            if bestSwitch[(ep, eh2, l_top, l_bottom)] == O[(ep, eh1, l_top, l_bottom)] and \
                            oBest[(vp, hChild1, l_top, l_bottom)]!=[(None, None, None, None)]:
                                bestSwitchLocations[(vp, hChild2, l_top, l_bottom)].extend\
                                (oBest[(vp, hChild1, l_top, l_bottom)])

    for key in bestSwitchLocations:
        if bestSwitchLocations[key][0] == (None, None, None, None):
            bestSwitchLocations[key] = bestSwitchLocations[key][1:]
    # Add the costs of each event to the corresponding eventsDict entry
    for key in eventsDict:
        eventsDict[key].append(Minimums[key])

    # Use findPath and findBestRoots to construct the DTLOR graph dictionary
    #this finds the optimal solutions (mappings)
    treeMin = findBestRoots(parasiteTree, Minimums)
    DTLOR = findPath(treeMin, eventsDict, {})
    for key in Score.keys():
        if not key in DTLOR:
            del Score[key]

    DTLOR, numRecon = addScores(treeMin, DTLOR, Score)

    return DTLOR, numRecon


def preorderDTLORsort(DTLOR, ParasiteRoot):
    """This takes in a DTL reconciliation graph and parasite root and returns 
    a sorted list, orderedKeysL, that is ordered by level from largest to 
    smallest, where level 0 is the root and the highest level has tips."""

    keysL = Greedy.orderDTL(DTLOR, ParasiteRoot)
    uniqueKeysL = Greedy.sortHelper(DTLOR, keysL)
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
            #TODO
            vertices = root[0] #format (vp, vh, l, l)
            if root[1] == 0:
                parentsDict[vertices] = ScoreDict[vertices]
            for n in range(len(DTLDict[vertices])-1): #the last item is always the cost
                _,child1,child2,oldScore = DTLDict[vertices][n]
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
    and returns a list of the minimum cost reconciliation tree roots"""
    treeTops = []
    for key in MinimumDict:
        #pick out the set of keys in C (MinimumDict) that correspond to mapping 
        #the root of parasite to any place in host tree
        #and any locus combinations
        if key[0] == Parasite['pTop'][1]:
            treeTops.append(key)
    treeMin = []
    for pair in treeTops:
        #pick out the set of keys that gives the minimum of C (cost of optimal solution)
        #note there might be multiple equally optimal mappings
        if MinimumDict[pair] == min([MinimumDict[root] for root in treeTops]):
            treeMin.append(pair)
    return treeMin

def findPath(tupleList, eventDict, uniqueDict):
    """Takes as input tupleList, a list of minimum reconciliation cost mappings 
     of the form (vp, vh, l_top, l_bottom), eventDict, the dictionary of events and 
     children for each node, and uniqueDict, the dictionary of unique vertex mappings. 
     This returns the completed DTL graph as a Dictionary"""
    for mapping in tupleList:
        if not mapping in uniqueDict:
            uniqueDict[mapping] = eventDict[mapping]
        for event in eventDict[mapping][:-1]:
            for location in event:
                if type(location) is tuple and location != (None, None, None, None):
                    findPath([location], eventDict, uniqueDict)
    return uniqueDict

def reconcile(fileName, D, T, L):
    """Takes as input a newick file, FileName, a dupliction cost, a transfer 
    cost, and a loss cost. This uses newickFormatReader to extract the host 
    tree, parasite tree and tip mapping from the file and then calls DP to 
    return the DTL reconciliation graph of the provided newick file"""
    host, paras, phi = newickFormatReader.getInput(fileName)
    return DP(host, paras, phi, D, T, L)
