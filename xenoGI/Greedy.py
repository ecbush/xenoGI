# Greedy.py
#Ivy Liu  
# Modified from code by
# Srinidhi Srinivasan, Juliet Forman
# June 2015

# This file contains the functions for finding all optimal reconciliations 
# using a DTL graph, as well as the scores for each of the events using the
# vertex-based DP algorithm. The main function in this file is called Greedy
# and the remaining functions are helper functions that are used by Greedy.

import copy 


def findRoot(Tree):
    """This function takes in a parasiteTree and returns a string with the 
    name of the root vertex of the tree"""

    if 'pTop' in Tree:
        return Tree['pTop'][1]
    return Tree['hTop'][1] 


def orderDTLOR(DTLOR, ParasiteRoot):
    """This function takes in a DTLOR graph and the ParasiteRoot. It outputs a 
    list, keysL, that contains tuples. Each tuple has two elements. The first
    is a mapping node of the form (p, h, l, l), where p is a parasite node and h is 
    a host node, l is the syntenic location of p. The second element is a level representing the depth of that
    mapping node within the tree."""

    keysL = []
    topNodes = []
    for key in DTLOR:
        if key[0] == ParasiteRoot:
            topNodes.append(key)
    for vertex in topNodes:
        keysL.extend(orderDTLOR_Roots(DTLOR, vertex, 0))
    return keysL


def orderDTLOR_Roots(DTLOR, vertex, level):
    """This function takes a DTL graph, one node, a vertex, of the DTL graph, 
    and level, and returns a list, keysL, that contains tuples. Each tuple has
    two elements. The first is a mapping node of the form (p, h, l, l), where p is a
    parasite node, h is a host node and l is the syntenic location of p. 
    The second element is a level representing the depth of that mapping node 
    within the tree. This function adds the input vertex to keysL and recurses 
    on its children."""

    keysL = []
    # Loop through each event associated with key in DTL
    for i in range(len(DTLOR[vertex]) - 1):
        event = DTLOR[vertex][i]
        child1 = event[1]
        child2 = event[2]
        keysL = keysL + [(vertex, level)]
        if child1[0] is not None:
            keysL.extend(orderDTLOR_Roots(DTLOR, child1, level + 1))
        if child2[0] is not None:
            keysL.extend(orderDTLOR_Roots(DTLOR, child2, level + 1)) 
    return keysL


def sortHelper(DTLOR, keysL):
    """This function takes in a list keysL,and a DTLOR graph, and returns a new 
    list, uniqueKeysL that has removed duplicate mapping nodes from keysL.
    This function chooses the highest level for each mapping node because we 
    are using the bottom-up approach."""
    
    uniqueKeysL = []
    for key in DTLOR:
        maxLevel = float("-inf")
        for element in keysL:
            if key == element[0]:
                if element[1] > maxLevel:
                    maxLevel = element[1]
        uniqueKeysL.append((key, maxLevel))
    return uniqueKeysL


def postorderDTLsort(DTLOR, ParasiteRoot):
    """This function takes in a DTL graph and ParasiteRoot, and returns a 
    sorted list, orderedKeysL, that is ordered by level from largest to 
    smallest, where level 0 is the root and the highest level are tips."""

    keysL = orderDTLOR(DTLOR, ParasiteRoot)
    uniqueKeysL = sortHelper(DTLOR, keysL)
    orderedKeysL = []
    levelCounter = 0
    while len(orderedKeysL) < len(uniqueKeysL):
        for mapping in uniqueKeysL:
            if mapping[-1] == levelCounter:
                orderedKeysL = [mapping] + orderedKeysL
        levelCounter += 1
    return orderedKeysL


def bookkeeping(DTLOR, ParasiteTree):
    """This function takes as inputs a DTLOR graph and ParasiteTree, and then 
    records what the max is at each mapping node and which event the max came 
    from. It outputs a dictionary BSFHMap by looping through the keys in a 
    sorted list of mapping nodes and finding the max score at each mapping 
    node and event node. BSFHMap has keys of the form (p, h, l, l) which are the 
    mapping nodes, and values which are lists where the first element is an 
    event with the max score, and the last element is the maxScore."""

    # Example: BSFHMap = {(mapping node): [['event', (p, h, l, l), (p, h, l, l), score], maxScore]}
    # Example: BSFHEvent = {(event node): max}

    BSFHMap = {}
    BSFHEvent = {}
    ParasiteRoot = findRoot(ParasiteTree)
    orderedKeysL = postorderDTLsort(DTLOR, ParasiteRoot)   
    for key in orderedKeysL:
        mapNode = key[0]
        # Check if the key is a tip:
        if DTLOR[mapNode][0][0] == 'C':              
            BSFHMap[mapNode] = [DTLOR[mapNode][0], DTLOR[mapNode][0][-1]]
        # If the key isn't a tip:
        else: 
            # Initialize counter
            maxScore = float("-inf")  
            # Initialize variable to keep track of where max came from
            maxEvent = [] 
            # Iterate through the events associated with the key node
            for i in range(len(DTLOR[mapNode]) - 1):  
                event = tuple(DTLOR[mapNode][i])
                BSFHEvent[event] = event[-1]
                if event[1] != (None, None, None, None):
                    BSFHEvent[event] = BSFHEvent[event]+BSFHMap[event[1]][-1]
                if event[2] != (None, None, None, None):
                    BSFHEvent[event] = BSFHEvent[event]+BSFHMap[event[2]][-1]
                # Check if current event has a higher score than current max
                if BSFHEvent[event] > maxScore:  
                    maxScore = BSFHEvent[event]  # Set new max score
                    maxEvent = list(event)   # Record where new max came from
            BSFHMap[mapNode] = [maxEvent, maxScore]  # Set BSFH value of key
    return BSFHMap


def TraceChildren(DTLOR, GreedyOnce, BSFHMap, key):
    """This function takes as input a DTLOR graph, a dictionary GreedyOnce, 
    containing the root of an optimal reconciliation, a BSFHMap dictionary,
    and a current key, and returns the optimal reconciliation and new DTLOR 
    graph, where the scores of the events in the reconciliation are reset to
    0."""

    resetDTLOR = {}  # The new DTL graph
    reset1DTLOR = {}  # This DTL graph deals with the recursive call on the child1
    reset2DTLOR = {}  # This DTL graph deals with the recursive call on the child2
    child1 = GreedyOnce[key][1]
    child2 = GreedyOnce[key][2]
    if child1 != (None, None, None, None):
        GreedyOnce[child1] = BSFHMap[child1][0][0:3]  # Add event to greedyOnce, without the score
        # This loop resets all the scores of events that have been used to 0
        for i in range(len(DTLOR[child1]) - 1):       
            if DTLOR[child1][i] == BSFHMap[child1][0]:
                newValue = DTLOR[child1]
                newValue[i][-1] = 0
                reset1DTLOR[child1] = newValue
        # This recursive call updates GreedyOnce and the DTL graph
        newGreedyOnce, DTLOR1 = TraceChildren(DTLOR, GreedyOnce, BSFHMap, child1) 
        reset1DTLOR.update(DTLOR1)
        GreedyOnce.update(newGreedyOnce)
    if child2 != (None, None, None, None):
        GreedyOnce[child2] = BSFHMap[child2][0][0:3]
        for i in range(len(DTLOR[child2]) - 1): 
            if DTLOR[child2][i] == BSFHMap[child2][0]:
                newValue = DTLOR[child2]
                newValue[i][-1] = 0
                reset2DTLOR[child2] = newValue      
        newGreedyOnce, DTLOR2 = TraceChildren(DTLOR, GreedyOnce, BSFHMap, child2)
        reset2DTLOR.update(DTLOR2)
        GreedyOnce.update(newGreedyOnce)
    resetDTLOR.update(reset1DTLOR)
    resetDTLOR.update(reset2DTLOR)
    return GreedyOnce, resetDTLOR


def greedyOnce(DTLOR, ParasiteTree):
    """This function takes the DTL graph and the ParasiteTree as inputs, and 
    calls bookkeeping to find the dictionary BSFHMap. It returns the 
    reconciliation with the highest score in a dictionary called GreedyOnce, 
    and also resets to 0 the scores of the mapping nodes in the optimal 
    reconciliation. The returned dictionary will have keys which are the 
    mapping nodes in the optimal reconciliation, and values of the form 
    (event, child1, child2)."""

    BSFHMap = bookkeeping(DTLOR, ParasiteTree)
    ParasiteRoot = findRoot(ParasiteTree)
    GreedyOnce = {}            # Initialize dictionary we will return
    bestKey = ()       # Variable to hold the key with the highers BSFH value
    bestScore = float("-inf")  # Variable to hold the highest BSFH value so far
    # Iterate trough all the keys (vertices) in BSFHMap
    for key in BSFHMap:   
        # Check if key has a score higher than bestScore and includes PRoot
        if BSFHMap[key][-1] > bestScore and key[0] == ParasiteRoot: 
            bestKey = key
            bestScore = BSFHMap[key][-1]
    # Set value in GreedyOnce of the best key we found
    GreedyOnce[bestKey] = BSFHMap[bestKey][0][0:3]                  
    
    # Reset score of the mapping node we used to 0:

    # Loop through the events associated with DTL
    for i in range(len(DTLOR[bestKey]) - 1):  
        # Check if the event matches the event that gave the best score
        if DTLOR[bestKey][i] == BSFHMap[bestKey][0]:           
            newEvent = DTLOR[bestKey][i][:-1] + [0]
            newValue = DTLOR[bestKey][:i] + [newEvent] + DTLOR[bestKey][i + 1:]
            DTLOR[bestKey] = newValue  # Set the score to 0
    newGreedyOnce, resetDTLOR = TraceChildren(DTLOR, GreedyOnce, BSFHMap, bestKey)
    GreedyOnce.update(newGreedyOnce)
    DTLOR.update(resetDTLOR)
    return GreedyOnce, DTLOR, bestScore


def Greedy(DTLOR, ParasiteTree):
    """This function takes as input a DTL graph and a ParasiteTree, and 
    returns TreeList, a list of dictionaries, each of which represent one of 
    the optimal reconciliations. This function runs till all the scores have 
    been collected from the DTL graph."""
    scores = []  # List of reconciliation scores
    currentDTLOR = copy.deepcopy(DTLOR)
    rec = [] # List of reconciliations
    collected = True
    while collected:
        # Call greedyOnce if all the points have not been collected yet
        oneTree, currentDTLOR, score = greedyOnce(currentDTLOR, ParasiteTree)
        scores.append(score) 
        rec.append(oneTree)
        collected = False  # Set not0 to False
        # Iterate to see if more points need to be collected
        for key in currentDTLOR:
            for i in range(len(currentDTLOR[key])-1):
                if currentDTLOR[key][i][-1] != 0:
                    collected = True
    return scores, rec
