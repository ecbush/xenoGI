import numpy,glob,struct

class Score:

    def __init__(self):
        '''Create an object for storing scores.'''
        
        self.endNodesToEdgeD = {}
        self.scoreD = {}

    def initializeDataAttributes(self,blastFilePath,geneNames):
        '''This method takes a new, empty object and fills the data attributes
by reading through blast files to identify pairs of genes with
edges. Also creates the arrays for storing scores, initialized to 0.'''
        self.fillEndNodesToEdgeD(blastFilePath,geneNames)
        self.numEdges=len(self.endNodesToEdgeD)
        self.initializeScoreD()

    def fillEndNodesToEdgeD(self,blastFilePath,geneNames):
        '''Run through blast files, finding all pairs of genes with signicant similarity. Use these to fill endNodesToEdgeD.'''

        blastFnL=glob.glob(blastFilePath)
        
        # Run through blast files getting list of all genes to compare. The
        # blast files have redundancies (ie , g1,g2 and also g2,g1). Keep
        # comparisons to do in set, and don't add if we already have it in
        # either order.
        edgeNum=0
        for fn in blastFnL:
            f = open(fn,'r')
            while True:
                s = f.readline()
                if s=='':
                    break
                L = s.split('\t')
                if len(L) != 12: # we only want lines with 12 columns
                    continue

                g1 = geneNames.nameToNum(L[0])
                g2 = geneNames.nameToNum(L[1])

                # make sure g1 is the lower gene number (we do this so
                # we only store in one orientation in dict).
                if g1 > g2: g2,g1 = g1,g2

                if (g1,g2) not in self.endNodesToEdgeD:
                    # we haven't got it yet, add it
                    self.endNodesToEdgeD[(g1,g2)] = edgeNum
                    edgeNum += 1
            f.close()

    def initializeScoreD(self):
        '''Create the arrays for storing scores.'''
        
        self.scoreD['rawSc'] = numpy.zeros(self.numEdges,dtype=numpy.float64)
        self.scoreD['normSc'] = numpy.zeros(self.numEdges,dtype=numpy.float64)
        self.scoreD['synSc'] = numpy.zeros(self.numEdges,dtype=numpy.float64)
        self.scoreD['coreSynSc'] = numpy.zeros(self.numEdges,dtype=numpy.float64)

    def isEdgePresentByEndNodes(self,g1,g2):
        '''See if g1 and g2 have an edge in our data structure. Return
boolean.'''
        if g1 > g2: g2,g1 = g1,g2
        if (g1,g2) in self.endNodesToEdgeD: return True
        else: return False
        
    def endNodesToEdge(self,g1,g2):
        '''Given two genes, return the number of the edge between them. If
there isn't any, return None.'''
        edge = self.endNodesToEdgeD[(g1,g2)]
        return edge

    def addScoreByEdge(self,edge,sc,scoreType):
        '''Given a score on an edge, store it in the score array corresponding
to scoreType.
        '''
        self.scoreD[scoreType][edge] = sc
        
    def addScoreByEndNodes(self,g1,g2,sc,scoreType):
        '''Given a score between two genes store it in the score array
corresponding to scoreType.'''

        if g1 > g2: g2,g1 = g1,g2 # make sure g1 is lower gene num

        if self.isEdgePresentByEndNodes(g1,g2):
            # we can only add scores for those edges that already exist
            edge = self.endNodesToEdge(g1,g2)
            self.addScoreByEdge(edge,sc,scoreType)
        # if the edge isn't already present, do nothing.
            
    def getScoreByEdge(self,edge,scoreType):
        '''Given an edge get the score corresponding to scoreType.'''
        return self.scoreD[scoreType][edge]

    def getScoreByEndNodes(self,g1,g2,scoreType):
        '''Given two genes get the score corresponding to scoreType.'''

        if g1 > g2: g2,g1 = g1,g2 # make sure g1 is lower gene num
        edge = self.endNodesToEdge(g1,g2)
        return self.getScoreByEdge(edge,scoreType)

    def iterateEdges(self):
        '''Returns an iterator which goes over edges by edge number.'''
        return range(self.numEdges)

    def iterateEdgesByEndNodes(self):
        '''Returns an iterator which goes over edges by end nodes.'''
        return self.endNodesToEdgeD.keys()

    def createNodeConnectL(self,geneNames):
        '''Create an attribute nodeConnectL. Index in this list corresponds to
gene. Value at that index is a list of the genes which a given gene
connects to. This attribute is not saved in our file formats. It must
be recalculated before it will be used (e.g. in family formation).'''

        self.nodeConnectL = [[] for gn in geneNames.nums]

        # loop over endNodesToEdgeD populating nodeConnectL
        for gn1,gn2 in self.endNodesToEdgeD:
            self.nodeConnectL[gn1].append(gn2)
            self.nodeConnectL[gn2].append(gn1)

    def getConnectionsGene(self,gene):
        '''Return a list containing all the genes connected to gene.'''
        return self.nodeConnectL[gene]

    def createNodeEdgeL(self,geneNames):
        '''Create an attribute nodeEdgeL. Index in this list corresponds to
gene. Value at that index is a list of the edges which a give gene
connects to. This attribute is not saved in our file formats. It must
be recalculated before it will be used.'''

        self.nodeEdgeL = [[] for gn in geneNames.nums]

        # loop over endNodesToEdgeD populating nodeEdgeL
        for gn1,gn2 in self.endNodesToEdgeD:
            edge = self.endNodesToEdge(gn1,gn2)
            self.nodeEdgeL[gn1].append(edge)
            self.nodeEdgeL[gn2].append(edge)

    def getConnectionsEdge(self,gene):
        '''Return a list containing all the edges connected to gene.'''
        return self.nodeEdgeL[gene]

    def createEdgeToEndNodeL(self):
        '''Create an attribute edgeToEndNodeL where the index is edge number,
and the value is a tuple giving the two genes on either end of the
edge. This attribute is not saved in our file formats. It must be
recalculated before it will be used.
        '''
        self.edgeToEndNodeL = [None for x in range(len(self.endNodesToEdgeD))]
        
        for endPairT in self.endNodesToEdgeD:
            edge = self.endNodesToEdgeD[endPairT]
            self.edgeToEndNodeL[edge] = endPairT

    def getEndNodesByEdge(self,edge):
        '''Given and edge, return the numbers for the two genes on either end.'''
        return self.edgeToEndNodeL[edge]
        
    def __eq__(self,other):
        '''Compare two Score objects to see if they have the same values.'''
        
        c1 = self.endNodesToEdgeD == other.endNodesToEdgeD
        c2 = all(self.scoreD['rawSc'] == other.scoreD['rawSc'])
        c3 = all(self.scoreD['normSc'] == other.scoreD['normSc'])
        c4 = all(self.scoreD['synSc'] == other.scoreD['synSc'])
        c5 = all(self.scoreD['coreSynSc'] == other.scoreD['coreSynSc'])

        if c1 and c2 and c3 and c4 and c5:
            return True
        else:
            return False
        
    def writeScoresText(self,geneNames,scoresFN):
        '''Save all edges (pairs of genes) to a tab delimited text file with a
    header line.
        '''

        scoreTypeL = ['rawSc','normSc','synSc','coreSynSc']

        # open file
        f=open(scoresFN,'w')

        # write header
        f.write("\t".join(['gene1','gene2','edge']+scoreTypeL)+'\n')

        for gene1Num,gene2Num in self.endNodesToEdgeD:
            edge = self.endNodesToEdgeD[(gene1Num,gene2Num)]
            outStrL=[]
            outStrL.append(geneNames.numToName(gene1Num))
            outStrL.append(geneNames.numToName(gene2Num))
            outStrL.append(str(edge))
            for scoreType in scoreTypeL:
                outStrL.append(format(self.getScoreByEndNodes(gene1Num,gene2Num,scoreType),".6f"))

            f.write("\t".join(outStrL)+'\n')

        f.close()

    def writeScoresBinary(self,scoresFN):
        '''Write scores to scoresFN as a binary file.'''

        scoreTypeL = ['rawSc','normSc','synSc','coreSynSc']
        f=open(scoresFN,'wb')

        # first 8 bytes are int, telling us how many edges
        f.write(self.numEdges.to_bytes(8,'little'))

        # write a block of bytes for each edge.
        # g1 g2 edge
        for g1,g2 in self.endNodesToEdgeD:
            edge = self.endNodesToEdgeD[(g1,g2)]
            f.write(g1.to_bytes(8,'little'))
            f.write(g2.to_bytes(8,'little'))
            f.write(edge.to_bytes(8,'little'))
            for scoreType in scoreTypeL:
                val = float(self.getScoreByEndNodes(g1,g2,scoreType))
                f.write(struct.pack('<d',val))
        
        f.close()
        
    def readScoresText(scoresFN,geneNames):
        '''Read scores from a text file of scores and use to create Score
        object
        '''

        scoresO = Score()
        
        # first pass over file, fill endNodesToEdgeD
        f = open(scoresFN,'r')
        f.readline() # skip header
        while True:
            s = f.readline()
            if s == '':
                break
            lineL=s.split('\t')
            g1 = geneNames.nameToNum(lineL[0])
            g2 = geneNames.nameToNum(lineL[1])
            edge = int(lineL[2])
            scoresO.endNodesToEdgeD[(g1,g2)] = edge
        f.close()
        
        # Set up the arrays
        scoresO.numEdges=len(scoresO.endNodesToEdgeD)
        scoresO.initializeScoreD()


        # 2nd pass, fill arrays
        f = open(scoresFN,'r')
        headerL=f.readline().split()
        scoreTypeL=headerL[3:] # first 3 cols are gene names + edge

        while True:
            s = f.readline()
            if s == '':
                break
            lineL=s.split('\t')
            edge = int(lineL[2])
            scoreL=lineL[3:]
            for i,sc in enumerate(scoreL):
                scoresO.scoreD[scoreTypeL[i]][edge] = float(sc)
        f.close()
        return scoresO

    def readScoresBinary(scoresFN):
        '''Read scores from a binary files.'''

        scoresO=Score()
        f=open(scoresFN,'rb')

        # read first 8 bytes to get numEdges
        b = f.read(8)
        scoresO.numEdges = int.from_bytes(b,'little')

        # initialize score arrays
        scoresO.initializeScoreD()
        
        # loop over blocks of bytes where each block is an edge
        for i in range(scoresO.numEdges):
            g1 = int.from_bytes(f.read(8),'little')
            g2 = int.from_bytes(f.read(8),'little')
            edge = int.from_bytes(f.read(8),'little')
            
            scoresO.endNodesToEdgeD[(g1,g2)] = edge

            # fill array
            scoresO.scoreD['rawSc'][edge] = struct.unpack('<d',f.read(8))[0]
            scoresO.scoreD['normSc'][edge] = struct.unpack('<d',f.read(8))[0]
            scoresO.scoreD['synSc'][edge] = struct.unpack('<d',f.read(8))[0]
            scoresO.scoreD['coreSynSc'][edge] = struct.unpack('<d',f.read(8))[0]            
            
        f.close()
        return scoresO
