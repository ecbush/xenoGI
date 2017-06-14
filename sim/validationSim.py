import sys,os
sys.path.append(os.path.join(sys.path[0],'../'))
import parameters,trees

def loadLog(logFN):
    '''Load a simulation log, and store the events in a tuple. index is
node. value is all the events that happened on that branch. We assume
the events at each node come in order in the file.'''

    # load it all into a dict keyed by node
    nodeEventD={}
    f=open(logFN,'r')
    while True:
        s=f.readline()
        if s=='':
            break
        L=s.rstrip().split('\t')
        node=L[0]
        evType = L[1]
        geneT = tuple(L[2:])

        if node in nodeEventD:
            nodeEventD[node].append((evType,geneT))
        else:
            nodeEventD[node] = [(evType,geneT)]

    # convert over to a tuple, where index is node number

    

            
def createSubLeafT(tree):
    '''Create a tuple, where index corresponds to tree node. The element
at that index is a list of the leaves descended from that node.'''
    subtreeL=trees.createSubtreeL(tree)
    subtreeL.sort()
    subLeafL=[]
    for subtree in subtreeL:
        subLeafL.append(tuple(trees.leafList(subtree)))
    subLeafT = tuple(subLeafL)
    return subLeafT
    
    
if __name__ == "__main__":

    paramFN=sys.argv[1]
    paramD = parameters.loadParametersD(paramFN)

    tree,strainStr2NumD,strainNum2StrD = trees.readTree(paramD['treeFN'])

    subLeafT = createSubLeafT(tree)

    
