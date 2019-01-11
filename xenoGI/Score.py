import numpy,glob,os,struct,statistics

class Score:

    def __init__(self):
        '''Create an object for storing scores.'''
        
        self.endNodesToEdgeD = {}
        self.numEdges = 0 # to start out
        self.strainPairScoreLocationD = {}
        self.scoreD = {}

    def initializeDataAttributes(self,blastFilePath,geneNames,strainStr2NumD):
        '''This method takes a new, empty object and fills the data attributes
by reading through blast files to identify pairs of genes with
edges. Also creates the arrays for storing scores, initialized to 0.'''
        self.fillEndNodesToEdgeD(blastFilePath,geneNames,strainStr2NumD)
        self.numEdges=len(self.endNodesToEdgeD)
        self.initializeScoreD()

    def fillEndNodesToEdgeD(self,blastFilePath,geneNames,strainStr2NumD):
        '''Run through blast files, finding all pairs of genes with signicant
similarity. Use these to fill endNodesToEdgeD. Also keep track of the
edge numbers associated with particular strain pairs and save in
strainPairScoreLocationD.
        '''

        # get blast files organized by strain pairs
        
        blastFnByPairD = self.getBlastFnByPairD(blastFilePath,strainStr2NumD)

        edgeNum=0
        for strainPair in blastFnByPairD:
            pairT = blastFnByPairD[strainPair]
            
            strainPairSt = edgeNum # scores for this strain pair will
                                   # begin at edgeNum

            # pairT may have two blast files (e.g for strain 1 vs 2
            # and also for 2 vs 1), or may only have one (e.g. for
            # strain 1 vs strain 1).
            for fn in pairT:
            
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

                    # because we have blast files going both ways
                    # (e.g. strain 1 vs strain 2 and also strain 2 vs
                    # strain 1), we check if we've added this
                    # particular gene pair already. If so, no need to
                    # do so again.
                    
                    # make sure g1 is the lower gene number (we do this so
                    # we only store in one orientation in dict).
                    if g1 > g2: g2,g1 = g1,g2

                    if (g1,g2) not in self.endNodesToEdgeD:
                        # we haven't got it yet, add it
                        self.endNodesToEdgeD[(g1,g2)] = edgeNum
                        edgeNum += 1
                f.close()
                
            strainPairEnd = edgeNum

            # Store the range of edge numbers for this
            # strainPair. Keys for this dict are tuples of train
            # number, e.g. (1,2) Always lower number first.
            self.strainPairScoreLocationD[strainPair] = (strainPairSt,strainPairEnd)
            

    def getBlastFnByPairD(self,blastFilePath,strainStr2NumD):
        '''Get the set of blast files and organize by the pair of strains
compared. Returns a dict keyed by tuple of strain number
(e.g. (1,2)). The value is a list with all the blast files comparing
two strains. For example one list might contain
['blast/E_albertii-E_coli_K12.out',
'blast/E_coli_K12-E_albertii.out']. For cases where a strain is
compared against itself then we have only one file name in the list.
        '''
        blastFnL=glob.glob(blastFilePath)

        blastFnByPairD ={}
        
        for fileStr in blastFnL:
            strain1,strain2 = os.path.splitext(os.path.split(fileStr)[-1])[0].split('-')
            strainNum1 = strainStr2NumD[strain1]
            strainNum2 = strainStr2NumD[strain2]            
            key = tuple(sorted([strainNum1,strainNum2]))
            if key in blastFnByPairD:
                # we've already looked at the one that had them in the other order
                blastFnByPairD[key].append(fileStr)
            else:
                blastFnByPairD[key] = [fileStr]
            
        return blastFnByPairD
                
    def initializeScoreD(self):
        '''Create the arrays for storing scores.'''
        
        self.scoreD['rawSc'] = numpy.zeros(self.numEdges,dtype=numpy.float64)
        self.scoreD['synSc'] = numpy.zeros(self.numEdges,dtype=numpy.float64)
        self.scoreD['coreSynSc'] = numpy.zeros(self.numEdges,dtype=numpy.float64)

    def getStrainPairs(self):
        '''Return all the strain pairs associated with this scores object as a
list of tuples. Each tuple will be two strains in numeric form.'''
        return sorted(self.strainPairScoreLocationD.keys())
        
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

    def iterateScoreByStrainPair(self,strainPair,scoreType):
        '''Returns an iterator which gives all scores of type scoreType associated with a
particular strainPair.'''
        stInd,endInd = self.strainPairScoreLocationD[strainPair]
        for edge in range(stInd,endInd):
            yield self.getScoreByEdge(edge,scoreType)

    def __eq__(self,other):
        '''Determine if other has the same values as self. Useful mainly for debugging.'''

        if self.numEdges != other.numEdges:
            return False

        if sorted(self.endNodesToEdgeD.items()) != sorted(other.endNodesToEdgeD.items()):
            return False

        if sorted(self.strainPairScoreLocationD.items()) != sorted(other.strainPairScoreLocationD.items()):
            return False

        if sorted(self.scoreD.keys()) != sorted(other.scoreD.keys()):
            return False
        
        for scoreType in self.scoreD.keys():
            if not all(self.scoreD[scoreType] == other.scoreD[scoreType]):
                return False

        # we made it, they're the same
        return True
        
    def writeScoresText(self,geneNames,scoresFN):
        '''Save all edges (pairs of genes) to a tab delimited text file with a
    header line.
        '''

        scoreTypeL = ['rawSc','synSc','coreSynSc']

        # open file
        f=open(scoresFN,'w')

        # Specify the index ranges for different strain pairs
        f.write("# Index range of scores for different strain pairs: strain1 strain2 stInd endInd.\n")
        for key in self.strainPairScoreLocationD:
            strain1,strain2 = key
            stInd,endInd = self.strainPairScoreLocationD[key]
            f.write("\t".join(map(str,[strain1,strain2,stInd,endInd]))+"\n")
            
        # write header
        f.write("# Scores: "+"\t".join(['gene1','gene2','edge']+scoreTypeL)+'\n')

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

    def readScoresText(scoresFN,geneNames):
        '''Read scores from a text file of scores and use to create Score
        object
        '''

        scoresO = Score()
        f = open(scoresFN,'r')

        # load the index ranges for different strain pairs
        f.readline() # skip first header/separator line
        
        maxEndInd = 0 # the biggest endInd we've seen so far
        while True:
            s = f.readline()
            if s[0] == '#':
                # we've hit the separator. scores begin on next line
                break
            lineL=s.split('\t')
            strain1 = int(lineL[0])
            strain2 = int(lineL[1])
            stInd = int(lineL[2])
            endInd = int(lineL[3])

            if endInd > maxEndInd:
                maxEndInd = endInd
            
            scoresO.strainPairScoreLocationD[(strain1,strain2)] = (stInd,endInd)

        # now we move to loading scores. The number of edges will be maxEndInd

        # Set up the arrays
        scoresO.numEdges = maxEndInd
        scoresO.initializeScoreD()

        separaterLineL = s.split()
        scoreTypeL=separaterLineL[5:] # the types of scores are listed in this separator line
        
        while True:
            s = f.readline()
            if s == '':
                break
            lineL=s.split('\t')

            g1 = geneNames.nameToNum(lineL[0])
            g2 = geneNames.nameToNum(lineL[1])
            edge = int(lineL[2])
            scoresO.endNodesToEdgeD[(g1,g2)] = edge

            scoreL=lineL[3:]
            for i,sc in enumerate(scoreL):
                scoresO.addScoreByEdge(edge,float(sc),scoreTypeL[i])

        f.close()
        return scoresO

    def writeScoresBinary(self,scoresFN):
        '''Write scores to scoresFN as a binary file.'''

        scoreTypeL = ['rawSc','synSc','coreSynSc']
        f=open(scoresFN,'wb')

        # first 8 bytes are int, telling us how many edges
        f.write(self.numEdges.to_bytes(8,'little'))

        # next 8 bytes are int telling us how many strain pairs there are
        numStrainPairs = len(self.strainPairScoreLocationD)
        f.write(numStrainPairs.to_bytes(8,'little'))

        # write a block of bytes for each strain pair location
        for key in self.strainPairScoreLocationD:
            strain1,strain2 = key
            stInd,endInd = self.strainPairScoreLocationD[key]
            f.write(strain1.to_bytes(8,'little'))
            f.write(strain2.to_bytes(8,'little'))
            f.write(stInd.to_bytes(8,'little'))
            f.write(endInd.to_bytes(8,'little'))
            
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
        
    def readScoresBinary(scoresFN):
        '''Read scores from a binary files.'''

        scoresO=Score()
        f=open(scoresFN,'rb')

        # read first 8 bytes to get numEdges
        b = f.read(8)
        scoresO.numEdges = int.from_bytes(b,'little')

        # read next 8 bytes to get the number of strain pairs
        b = f.read(8)
        numStrainPairs = int.from_bytes(b,'little')

        # read in the index locations where scores from different
        # strain pairs will be found
        for i in range(numStrainPairs):
            strain1 = int.from_bytes(f.read(8),'little')
            strain2 = int.from_bytes(f.read(8),'little')
            stInd = int.from_bytes(f.read(8),'little')
            endInd = int.from_bytes(f.read(8),'little')
            scoresO.strainPairScoreLocationD[(strain1,strain2)] = (stInd,endInd)
            
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
            scoresO.scoreD['synSc'][edge] = struct.unpack('<d',f.read(8))[0]
            scoresO.scoreD['coreSynSc'][edge] = struct.unpack('<d',f.read(8))[0]            
            
        f.close()
        return scoresO
    
    ## Below are functions that create various optional
    ## attributes. These attributes are not saved in the file format,
    ## and so must be calculated each time we intend to use them.

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

    def createAabrhScoreSummaryD(self,strainNumsL,aabrhL,geneNames):
        '''Given raw scores and set of all around best reciprocal hits,
    calculates the mean and standard deviation of scores and stores in a
    dictionary.'''

        strainNumsL.sort()
        
        # create dictionary keyed by species pair, (representing an upper
        # triangular matrix)
        spScoreD={}
        for i in range(len(strainNumsL)-1):
            strain1 = strainNumsL[i]
            for j in range(i+1,len(strainNumsL)):
                strain2 = strainNumsL[j]
                spScoreD[(strain1,strain2)]=[]

        # loop through aabrhL and populate
        for orthoT in aabrhL:
            spScoreD = self.addPairwiseScores(spScoreD,orthoT,geneNames)

        # get mean and standard deviation
        self.scoreSummaryD = {}
        for sp1,sp2 in spScoreD:
            mean = statistics.mean(spScoreD[(sp1,sp2)])
            std = statistics.stdev(spScoreD[(sp1,sp2)])
            self.scoreSummaryD[(sp1,sp2)] = (mean,std)
            self.scoreSummaryD[(sp2,sp1)] = (mean,std)

    def addPairwiseScores(self,spScoreD,orthoT,geneNames):
        '''Given a dictionary for storing pairwise scores, and ortholog set in
    orthoT, and a network of scores, scoresO, pull out all species pairs, and
    add score for each in appropriate place in spScoreD.'''

        for i in range(len(orthoT)-1):
            gene1 = orthoT[i]
            geneNum1=geneNames.nameToNum(gene1)
            sp1 = geneNames.nameToStrainNum(gene1)
            for j in range(i+1,len(orthoT)):
                gene2 = orthoT[j]
                geneNum2=geneNames.nameToNum(gene2)
                sp2 = geneNames.nameToStrainNum(gene2)
                sc = self.getScoreByEndNodes(geneNum1,geneNum2,'rawSc')
                key = tuple(sorted([sp1,sp2]))
                spScoreD[key].append(sc)
        return spScoreD
