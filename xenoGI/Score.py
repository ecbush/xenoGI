import numpy,glob,os,struct,statistics,ctypes
from multiprocessing import RawArray
from . import blast

class Score:

    def __init__(self):
        '''Create an object for storing scores.'''
        
        self.endNodesToEdgeD = {}
        self.numEdges = 0 # to start out
        self.strainPairScoreLocationD = {}
        self.scoreD = {}

    def initializeDataAttributes(self,blastFnL,paramD,strainNamesT,genesO):
        '''This method takes a new, empty object and fills the data attributes
by reading through blast files to identify pairs of genes with
edges. Also creates the array for storing raw scores, initialized to
0. Arrays for other score types will be created later. The genes
included here depend on what is found in the blast files in blastFnL.

        '''
        self.fillEndNodesToEdgeD(blastFnL,paramD,strainNamesT,genesO)
        self.numEdges=len(self.endNodesToEdgeD)
        self.initializeScoreArray('rawSc')

    def fillEndNodesToEdgeD(self,blastFnL,paramD,strainNamesT,genesO):
        '''Run through blast files, finding all pairs of genes with signicant
similarity. Use these to fill endNodesToEdgeD. Also keep track of the
edge numbers associated with particular strain pairs and save in
strainPairScoreLocationD.
        '''

        # get blast files organized by strain pairs
        evalueThresh = paramD['evalueThresh']
        alignCoverThresh = paramD['alignCoverThresh']
        percIdentThresh =  paramD['percIdentThresh']
        blastFileJoinStr = paramD['blastFileJoinStr']
        
        blastFnByPairD = self.getBlastFnByPairD(blastFnL,blastFileJoinStr,strainNamesT)

        # get edges
        edgeNum=0
        for g in genesO.iterGenes(): # first every gene vs itself
            if genesO.numToStrainName(g) in strainNamesT:
                # but only do it for genes that occur in a strain in strainNamesT
                self.endNodesToEdgeD[(g,g)] = edgeNum
                edgeNum += 1
            
        # now get edges out of blast files
        for strainPair in blastFnByPairD:
            pairT = blastFnByPairD[strainPair]
            
            strainPairSt = edgeNum # scores for this strain pair will
                                   # begin at edgeNum

            # pairT may have two blast files (e.g for strain 1 vs 2
            # and also for 2 vs 1), or may only have one (e.g. for
            # strain 1 vs strain 1).
            for fn in pairT:
                for g1,g2,evalue,alCov,pident,score,alLen in blast.parseBlastFile(fn,evalueThresh,alignCoverThresh,percIdentThresh):
                    
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
                    
                
            strainPairEnd = edgeNum

            # Store the range of edge numbers for this
            # strainPair. Keys for this dict are tuples of train
            # number, e.g. (1,2) Always lower number first.
            self.strainPairScoreLocationD[strainPair] = (strainPairSt,strainPairEnd)
            

    def getBlastFnByPairD(self,blastFnL,blastFileJoinStr,strainNamesT):
        '''Get the set of blast files and organize by the pair of strains
compared. Returns a dict keyed by tuple of strain number
(e.g. (1,2)). The value is a list with all the blast files comparing
two strains. For example one list might contain
['blast/E_albertii-E_coli_K12.out',
'blast/E_coli_K12-E_albertii.out']. For cases where a strain is
compared against itself then we have only one file name in the list.
        '''
        blastFnByPairD ={}
        
        for fileStr in blastFnL:
            strainName1,strainName2 = os.path.splitext(os.path.split(fileStr)[-1])[0].split(blastFileJoinStr)

            if strainName1 not in strainNamesT or strainName2 not in strainNamesT:
                # not a blast file we want, since it has someithing else besides our strains
                continue
            
            key = tuple(sorted([strainName1,strainName2]))
            if key in blastFnByPairD:
                # we've already looked at the one that had them in the other order
                blastFnByPairD[key].append(fileStr)
            else:
                blastFnByPairD[key] = [fileStr]
            
        return blastFnByPairD
                
    def initializeScoreArray(self,scoreType):
        '''Create array for storing scores.'''
        self.scoreD[scoreType] = numpy.zeros(self.numEdges,dtype=ctypes.c_double)

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

    def getStrains(self):
        '''Return a set of the strains used in this scores object (in number
form).'''
        strainS = set()
        for key in self.strainPairScoreLocationD:
            strain1,strain2 = key
            strainS.add(strain1)
            strainS.add(strain2)
        return strainS
        
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
        
    def writeScoresText(self,strainNamesT,scoreTypeL,scoresFN,genesO,geneInfoFN):
        '''Save all edges (pairs of genes) to a tab delimited text file with a
    header line.
        '''
        # open file
        f=open(scoresFN,'w')

        # Specify the index ranges for different strain pairs
        f.write("# Index range of scores for different strain pairs: strain1 strain2 stInd endInd.\n")
        for key in self.strainPairScoreLocationD:
            strainName1,strainName2 = key
            strainNum1 = strainNamesT.index(strainName1)
            strainNum2 = strainNamesT.index(strainName2)

            stInd,endInd = self.strainPairScoreLocationD[key]
            f.write("\t".join(map(str,[strainNum1,strainNum2,stInd,endInd]))+"\n")

        genesO.initializeGeneNumToNameD(geneInfoFN,strainNamesT)
            
        # write header
        f.write("# Scores: "+"\t".join(['gene1','gene2','edge']+scoreTypeL)+'\n')

        for gene1Num,gene2Num in self.endNodesToEdgeD:
            edge = self.endNodesToEdgeD[(gene1Num,gene2Num)]
            outStrL=[]
            outStrL.append(genesO.numToName(gene1Num))
            outStrL.append(genesO.numToName(gene2Num))
            outStrL.append(str(edge))
            for scoreType in scoreTypeL:
                outStrL.append(format(self.getScoreByEndNodes(gene1Num,gene2Num,scoreType),".6f"))

            f.write("\t".join(outStrL)+'\n')

        f.close()

    def readScoresText(strainNamesT,scoreTypeL,scoresFN):
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
            strainNum1 = int(lineL[0])
            strainNum2 = int(lineL[1])
            stInd = int(lineL[2])
            endInd = int(lineL[3])
            strainName1 = strainNamesT[strainNum1]
            strainName2 = strainNamesT[strainNum2]

            if endInd > maxEndInd:
                maxEndInd = endInd

            key = tuple(sorted([strainName1,strainName2]))
            scoresO.strainPairScoreLocationD[key] = (stInd,endInd)

        # now we move to loading scores. The number of edges will be maxEndInd

        # Set up the arrays
        scoresO.numEdges = maxEndInd
        for scoreType in scoreTypeL:
            scoresO.initializeScoreArray(scoreType)

        separaterLineL = s.split()
        scoreTypeL=separaterLineL[5:] # the types of scores are listed in this separator line
        
        while True:
            s = f.readline()
            if s == '':
                break
            lineL=s.split('\t')

            g1 = int(lineL[0].split('_')[0])
            g2 = int(lineL[1].split('_')[0])
            
            edge = int(lineL[2])
            scoresO.endNodesToEdgeD[(g1,g2)] = edge

            scoreL=lineL[3:]
            for i,sc in enumerate(scoreL):
                scoresO.addScoreByEdge(edge,float(sc),scoreTypeL[i])

        f.close()
        return scoresO

    def writeScoresBinary(self,strainNamesT,scoreTypeL,scoresFN):
        '''Write scores to scoresFN as a binary file.'''

        f=open(scoresFN,'wb')

        # first 8 bytes are int, telling us how many edges
        f.write(self.numEdges.to_bytes(8,'little'))

        # next 8 bytes are int telling us how many strain pairs there are
        numStrainPairs = len(self.strainPairScoreLocationD)
        f.write(numStrainPairs.to_bytes(8,'little'))

        # write a block of bytes for each strain pair location
        for key in self.strainPairScoreLocationD:
            strainName1,strainName2 = key
            strainNum1 = strainNamesT.index(strainName1)
            strainNum2 = strainNamesT.index(strainName2)
            stInd,endInd = self.strainPairScoreLocationD[key]
            f.write(strainNum1.to_bytes(8,'little'))
            f.write(strainNum2.to_bytes(8,'little'))
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
        
    def readScoresBinary(strainNamesT,scoreTypeL,scoresFN):
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
            strainNum1 = int.from_bytes(f.read(8),'little')
            strainNum2 = int.from_bytes(f.read(8),'little')
            stInd = int.from_bytes(f.read(8),'little')
            endInd = int.from_bytes(f.read(8),'little')
            strainName1 = strainNamesT[strainNum1]
            strainName2 = strainNamesT[strainNum2]
            key = tuple(sorted([strainName1,strainName2]))
            scoresO.strainPairScoreLocationD[key] = (stInd,endInd)
            
        # initialize score arrays
        for scoreType in scoreTypeL:
            scoresO.initializeScoreArray(scoreType)
        
        # loop over blocks of bytes where each block is an edge
        for i in range(scoresO.numEdges):
            g1 = int.from_bytes(f.read(8),'little')
            g2 = int.from_bytes(f.read(8),'little')
            edge = int.from_bytes(f.read(8),'little')
            
            scoresO.endNodesToEdgeD[(g1,g2)] = edge

            # fill array
            for scoreType in scoreTypeL:
                scoresO.scoreD[scoreType][edge] = struct.unpack('<d',f.read(8))[0]
            
        f.close()
        return scoresO
    
    ## Below are functions that create various optional
    ## attributes. These attributes are not saved in the file format,
    ## and so must be calculated each time we intend to use them.

    def createNodeConnectD(self):
        '''Create an attribute nodeConnectD. Keys corresponds to gene. Value
is a list of the genes which a given gene connects to. This attribute
is not saved in our file formats. It must be recalculated before it
will be used (e.g. in family formation).
        '''

        self.nodeConnectD = {}

        # loop over endNodesToEdgeD populating nodeConnectD
        for gn1,gn2 in self.endNodesToEdgeD:
            if gn1 not in self.nodeConnectD:
                self.nodeConnectD[gn1] = [gn2]
            else:
                self.nodeConnectD[gn1].append(gn2)
            if gn2 not in self.nodeConnectD:
                self.nodeConnectD[gn2] = [gn1]
            else:
                self.nodeConnectD[gn2].append(gn1)

    def getConnectionsGene(self,gene):
        '''Return a list containing all the genes connected to gene.'''
        if gene in self.nodeConnectD:
            return self.nodeConnectD[gene]
        else:
            return []
        
    def createNodeEdgeL(self,geneNamesO):
        '''Create an attribute nodeEdgeL. Index in this list corresponds to
gene. Value at that index is a list of the edges which a give gene
connects to. This attribute is not saved in our file formats. It must
be recalculated before it will be used.'''

        self.nodeEdgeL = [[] for gn in geneNamesO.iterGeneNums()]

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

    def createAabrhScoreSummaryD(self,strainNamesT,aabrhL,genesO):
        '''Given raw scores and set of all around best reciprocal hits,
    calculates the mean and standard deviation of scores and stores in a
    dictionary.'''

        # create dictionary keyed by species pair, (representing an upper
        # triangular matrix)
        spScoreD={}
        for i in range(len(strainNamesT)-1):
            strain1 = strainNamesT[i]
            for j in range(i+1,len(strainNamesT)):
                strain2 = strainNamesT[j]
                key = tuple(sorted([strain1,strain2]))
                spScoreD[key]=[]

        # loop through aabrhL and populate
        for orthoT in aabrhL:
            spScoreD = self.addPairwiseScores(spScoreD,orthoT,genesO)

        # get mean and standard deviation
        self.scoreSummaryD = {}
        for sp1,sp2 in spScoreD:
            mean = statistics.mean(spScoreD[(sp1,sp2)])
            std = statistics.stdev(spScoreD[(sp1,sp2)])
            self.scoreSummaryD[(sp1,sp2)] = (mean,std)
            self.scoreSummaryD[(sp2,sp1)] = (mean,std)

    def addPairwiseScores(self,spScoreD,orthoT,genesO):
        '''Given a dictionary for storing pairwise scores, and ortholog set in
    orthoT, and a network of scores, scoresO, pull out all species pairs, and
    add score for each in appropriate place in spScoreD.'''

        for i in range(len(orthoT)-1):
            gene1 = orthoT[i]
            sp1 = genesO.numToStrainName(gene1)
            for j in range(i+1,len(orthoT)):
                gene2 = orthoT[j]
                sp2 = genesO.numToStrainName(gene2)
                sc = self.getScoreByEndNodes(gene1,gene2,'rawSc')
                key = tuple(sorted([sp1,sp2]))
                spScoreD[key].append(sc)
        return spScoreD

    
class sharedScore:

    def __init__(self):
        '''Create an object for storing scores which will efficiently share
memory across processes.'''

        self.scoreD = {}
    
    def createArrays(self,scoresO,paramD):
        '''Creates a group of multiprocessing raw arrays for use in a
sharedScores object.
        '''

        hashArrayScaleFactor = paramD['hashArrayScaleFactor']
        
        numEdges = scoresO.numEdges
        self.hashArrayLen = hashArrayScaleFactor * numEdges
        
        ## Preliminary processing
        tempHashL=[[] for i in range(self.hashArrayLen)]
        maxGnNum=0 # get maximum numeric value used to represent a gene
        for gn1,gn2 in scoresO.iterateEdgesByEndNodes():
            edge = scoresO.endNodesToEdge(gn1,gn2)
            hs = self.hashByGenePair(gn1,gn2)
            tempHashL[hs].append((gn1,gn2,edge))
            
            if max(gn1,gn2) > maxGnNum:
                maxGnNum = max(gn1,gn2)
            
        # get length for collisionAr
        collisionArLen = 0
        maxEdgesPerHash = 0
        for L in tempHashL:
            L.sort() # sort on gn1
            if len(L) > 0:
                # add extra 1 because first block of this region will
                # hold the length of the region
                collisionArLen += len(L)
            if len(L) > maxEdgesPerHash:
                maxEdgesPerHash = len(L)
        self.collisionArLen = collisionArLen
                
        ## some checks
        
        # make sure data isn't too big. Assumes hash and gene arrays are c_uint32 (32
        # bit, unsigned)
        maxPossibleGeneOrEdgeNumber = 2**32-1
        if numEdges-1 > maxPossibleGeneOrEdgeNumber:
            raise ValueError("Error creating sharedScore object. Data set has too many score pairs for the data type we're using in our shared score arrays (c_uint32).")

        if maxEdgesPerHash > 2**16:
            raise ValueError("Error creating sharedScore object. The maximum number of score pairs per hash exceeds what lenOfRegionA can hold.")

        if maxGnNum > maxPossibleGeneOrEdgeNumber:
            raise ValueError("Error creating sharedScore object. At least one gene has been given a number too large for the array data type we are using.")


        ## create arrays

        rawScoreAr = RawArray(ctypes.c_double, numEdges)
        hasEdgeAr = RawArray(ctypes.c_bool, self.hashArrayLen) # tells if any edges at given hash
        hashAr = RawArray(ctypes.c_uint32, self.hashArrayLen) # stores index to col Ars
        lenOfRegionAr = RawArray(ctypes.c_uint16, self.hashArrayLen) # stores length of region holding stuff from a given hash
        colGn1Ar = RawArray(ctypes.c_uint32, self.collisionArLen)
        colGn2Ar = RawArray(ctypes.c_uint32, self.collisionArLen)
        colEdgeAr = RawArray(ctypes.c_uint32, self.collisionArLen)
        
        ## fill rawScoreAr
        for edge in scoresO.iterateEdges():
            sc = scoresO.getScoreByEdge(edge,'rawSc')
            rawScoreAr[edge] = sc

        ## fill hash and collision Ars
        colAr_i = 0 # index into col Ars, e.g. colGn1Ar
        for hs in range(len(tempHashL)):
            if tempHashL[hs] == []:
                hasEdgeAr[hs] = False
            else:
                hasEdgeAr[hs] = True
                lenOfRegionAr[hs] = len(tempHashL[hs])
                hashAr[hs] = colAr_i
                for gn1,gn2,edge in tempHashL[hs]:
                    colGn1Ar[colAr_i] = gn1
                    colGn2Ar[colAr_i] = gn2
                    colEdgeAr[colAr_i] = edge
                    colAr_i += 1

        self.insertArrays(rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,self.hashArrayLen)

    def hashByGenePair(self,gn1,gn2):
        '''Hash to an int, modulo array len.'''
        # python hash works better than cantor pairing func
        return hash((gn1,gn2)) % self.hashArrayLen

    def insertArrays(self,rawScoreAr,hasEdgeAr,hashAr,lenOfRegionAr,colGn1Ar,colGn2Ar,colEdgeAr,hashArrayLen):
        '''Attach the input arrays to self.'''
        self.scoreD['rawSc'] = rawScoreAr
        self.hasEdgeAr = hasEdgeAr
        self.hashAr = hashAr
        self.lenOfRegionAr = lenOfRegionAr
        self.colGn1Ar = colGn1Ar
        self.colGn2Ar = colGn2Ar
        self.colEdgeAr = colEdgeAr
        self.hashArrayLen = hashArrayLen

    def returnArrays(self):
        '''Return all our arrays.'''
        return self.scoreD['rawSc'],self.hasEdgeAr,self.hashAr,self.lenOfRegionAr,self.colGn1Ar,self.colGn2Ar,self.colEdgeAr,self.hashArrayLen
        
    def endNodesToEdge(self,gn1,gn2):
        '''Given two genes, return the number of the edge between them. If
there isn't any, return None.'''

        hs = self.hashByGenePair(gn1,gn2)

        if not self.hasEdgeAr[hs]:
            # no score pairs here
            return None
        else:
            # one or more gene pairs are stored in a region
            # corresponding to this hash. check if any of them are the
            # gene pair we're after.
            stInd = self.hashAr[hs]
            regionLen = self.lenOfRegionAr[hs]
            return self.searchCollisionArrays(gn1,gn2,stInd,regionLen)

    def searchCollisionArrays(self,gn1,gn2,stInd,regionLen):
        '''Search the collision arrays for a match to gn1,gn2. If found, return index,
    otherwise return None.
        '''

        # find the first occurence of gn1, if any
        # The colGn1Ar is sorted in numerical order between stInd and stInd+regionLen
        endInd = stInd + regionLen
        searchStart = self.binarySearchGn1Array(gn1,stInd,endInd)
        if searchStart == None:
            return None

        for i in range(searchStart,endInd):
            if self.colGn1Ar[i] != gn1:
                # it's no longer gn1, thus we won't find it
                return None
            elif self.colGn2Ar[i] == gn2:
                # self.colGn1Ar[i] == gn1 since we got this far
                return self.colEdgeAr[i]
        return None
        
    def binarySearchGn1Array(self,gn1,st,end):
        '''Find the first occurrence of gn1 in self.colGn1Ar using binary
search. (colGn1Ar can be assumed to be sorted over the region st to
end).
        '''
        while st < end:
            mid = (st + end) // 2
            if self.colGn1Ar[mid] < gn1:
                st = mid + 1
            elif self.colGn1Ar[mid] > gn1:
                end = mid
            elif mid > 0 and self.colGn1Ar[mid-1] == gn1:
                end = mid
            else:
                return mid

        return None

    def getScoreByEdge(self,edge,scoreType):
        '''Given an edge get the score corresponding to scoreType.'''
        return self.scoreD[scoreType][edge]

    def getScoreByEndNodes(self,g1,g2,scoreType):
        '''Given two genes get the score corresponding to scoreType.'''

        if g1 > g2: g2,g1 = g1,g2 # make sure g1 is lower gene num
        edge = self.endNodesToEdge(g1,g2)
        if edge == None:
            return None
        else:
            return self.getScoreByEdge(edge,scoreType)
