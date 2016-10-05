# functions for loading and manipulating phylogenetic trees

def middleComma(s):
    '''Find the 'middle' comma. This is the one which has had the same
number of ( as ) before it. Return the index.'''
    parenBalance=0
    for i in range(len(s)):
        if s[i]=='(':
            parenBalance+=1
        elif s[i]==')':
            parenBalance-=1
        elif s[i]==',' and parenBalance==0:
            return i
    return None

def newick2TupleTree(s):
    '''Given a Newick type string representation of a tree, return in
tuple tree format. If no branch length specified, makes it 0 length.'''
    if "(" not in s and ")" not in s:
        # its a tip
        name,brLen=s.split(":")
        brLen=float(brLen)
        return(name,(),(),brLen)
    else:
        # if there is a ':' on the right, outside the parens, split on
        # this for branch laength
        firstColonInd=len(s)-(s[::-1].index(":")+1)
        firstParenInd=len(s)-(s[::-1].index(")")+1)
        if firstColonInd>firstParenInd:
            brLen=float(s[firstColonInd+1:])
            s=s[:firstColonInd]
        else:
            brLen=0
        s=s[1:-1] # strip of outer perens
        mci=middleComma(s)
        leftTree=newick2TupleTree(s[:mci].strip()) # strip any leading or
        rightTree=newick2TupleTree(s[mci+1:].strip()) # trailing white space
        return('anc',leftTree,rightTree,brLen)

def nodeCount(tree):
    '''How many nodes in tree?'''
    if tree[1]==():
        return 1
    else:
        return 1+ nodeCount(tree[1]) + nodeCount(tree[2])

def leafList(tree):
    '''Return list of nodes in tree.'''
    if tree[1]==():
        return [tree[0]]
    else:
        return leafList(tree[1]) + leafList(tree[2])
    
    
def strTree2numTree(tree,counter):
    '''Given a tuple tree with nodes specified by strings, convert to
numbered nodes. Return tree with numbered nodes.'''
    if tree[1]==():
        return (counter,(),(),tree[3]),counter+1
    else:
        leftNumTree,counter=strTree2numTree(tree[1],counter)
        rightNumTree,counter=strTree2numTree(tree[2],counter)
        numTree=(counter,leftNumTree,rightNumTree,tree[3])
        return numTree,counter+1

def makeTreeD(tree1,tree2,treeD):
    '''Make a dictionary to convert from node names in tree1 to node names
in the identically shaped tree2.'''
    treeD[tree1[0]]=tree2[0]
    if tree1[1]==(): return
    else:
        makeTreeD(tree1[1],tree2[1],treeD)
        makeTreeD(tree1[2],tree2[2],treeD)
        return
    
    
def readTree(filename):
    '''Read Newick tree from file and convert it to tuple format.'''
    f=open(filename,"r")
    s=f.read().rstrip()
    f.close()
    if s[-1]==';': s=s[:-1] # strip off ';'
    stringTree=newick2TupleTree(s)
    counter=0
    numTree,counter=strTree2numTree(stringTree,counter)

    # make dictionaries for converting between number and string strain names
    strainStr2NumD={} # will be shorter due to collisions on 'anc'
    makeTreeD(stringTree,numTree,strainStr2NumD)
    strainNum2StrD={}
    makeTreeD(numTree,stringTree,strainNum2StrD)
    return numTree,strainStr2NumD,strainNum2StrD

def createSubtreeL(tree):
    '''Return a list containing all subtrees.'''
    if tree[1]==():
        return [tree]
    else:
        return [tree]+createSubtreeL(tree[1]) + createSubtreeL(tree[2])

def tupleTree2Newick(tree):
    '''Convert a four tuple based tree (root,left,right,branchLen) into a
newick formated string.'''
    if tree[1]==():
        return str(tree[0])+":"+str(tree[3])
    else:
        leftStr=tupleTree2Newick(tree[1])
        rightStr=tupleTree2Newick(tree[2])
        return "("+leftStr+","+rightStr+"):"+str(tree[3])

def writeTree(tree,fileName):
    '''Write tree to fileName (in newick format).'''
    f=open(fileName,"w")
    f.write("("+tupleTree2Newick(tree)+");")
    f.close()
