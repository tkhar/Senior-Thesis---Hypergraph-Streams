# Simplex Counting Experiments in Streamed Hypergraphs
# Testing out algorithms for Senior Thesis project

# Author: Haris Themistoklis
# Advisor: Amit Chakrabarti

# Date: 5/15/2021

import random
import copy
import statistics
import math

class HypergraphGenerator:
    def __init__(self, N, p):
        """
        This generator will be constructing random hypergraphs by selecting each 
        edge with probability p

        Args:
            N: the number of vertices
            p: the probability of having an edge
        """
        self.vertexNum = N
        self.hyperedgeProb = p

    def generate(self, filename):
        """
        Will generate the random hypergraph and save it in a file

        Args:
            filename: name of file to save produced hypergraph

        Returns:
            the number of edges in the constructed hypergraph
        """
        # Open the file we will write the hypergraph in
        file_writer = open(filename, "w")

        edgeCount = 0
        for i in range(1, self.vertexNum+1):
            for j in range(1, self.vertexNum+1):
                for k in range(1, self.vertexNum+1):
                    # ignore non-well ordered triples
                    if not(i < j and j < k):
                        continue
                    # keep hyperedge with probability p, simulating coin toss
                    tossCoin = random.random()
                    if tossCoin < self.hyperedgeProb:
                        file_writer.write(str(i)+" "+str(j)+" "+str(k)+"\n")
                        edgeCount = edgeCount +1
        
        print("Hypergraph generated! Edge Count = ", edgeCount)
        file_writer.close()
        
        return edgeCount


class Stream:
    def __init__(self, filename, N, M, MHS):
        """
        Initialize a stream of hyperedges

        Args:
            filename: the filename containing the hypergraph dataset
            N: the number of vertices in the hypergraph
            M: the number of hyperedges in the hypergraph
            MHS: the median hyperedge size in our hypergraph

        
        ASSUMPTIONS
            1. Each hyperedge is provided with all its elements being unique
            2. No duplicate hyperedges!
        """
        # file where the hyperedge dataset lives
        self.filename = filename              
        # getting a handle on this file                 
        self.file_reader = open(self.filename)
        # the number of vertices in this hypergraph
        self.n = N                       
        # the number of edges in this hypergraph
        self.m = M                     
        # the median hyperedge size  
        self.mhs = MHS                                        
        # the counts of the simplices we will be interested in approximately counting
        self.simplexCounts = [0 for i in range(self.mhs+3)]

    def done(self):
        self.file_reader.close()

    def CountAllSimplices(self,k=0):
        if k!=0:
            self.CountSimplices(k)
        else:
            for s in range(3, self.mhs+3):
                self.CountSimplices(s)
                self.file_reader.seek(0, 0)

    def CountSimplices (self, s):
        """
        Count the ground-truth number of s-simplices (simplices with s elements)
        We implement a naive algorithm that uses Theta(m) bits of space and O(m^2 n) time

        Args:
            s: the size of the simplex we want to count the occurences of
        """

        print("Counting simplices of size ",s)
        # s is outside of range of simplex we count
        if s <= 2 or s >= self.mhs + 3:
            return
        
        # Store all the edges in the hypergraph
        edgeSet = []
        edgeString = self.file_reader.readline()
        while edgeString:
            edgeSet.append(edgeString.split())
            edgeString = self.file_reader.readline()

        # This naive algorithm costs Theta(m) bits of memory
        self.memoryNaive = len(edgeSet)
        print("The naive method uses ", self.memoryNaive, " words of memory")

        progress = 0
        percentProgress = 0
        for e in edgeSet:
            # Keep track of progress
            progressed = False
            while progress >= percentProgress*len(edgeSet):
                percentProgress = percentProgress + 0.01
                progressed = True
            if progressed and (int(percentProgress*100))%20 == 0:
                print("Progress: ", int(percentProgress*100), "%")
            progress = progress + 1

            # Only look at edges with s-1 edges
            if len(e) != s-1:
                continue
            for v in range(1,self.n+1):
                if str(v) in e:
                    continue
                # Now we test if the proposed simplex is actually a simplex
                proposedSimplex = copy.deepcopy(e)
                proposedSimplex.append(str(v))

                # How many hyperedges exist in the proposed simplex?
                SubedgeCount = 0 
                for se in edgeSet:
                    if len(se) != s-1:
                        continue
                    if set(se).issubset(set(proposedSimplex)):
                        SubedgeCount = SubedgeCount + 1

                # Assuming that no duplicate hyperedges exist and that
                # no duplicate vertices can exist within a hyperedge, 
                # IF the proposed simplex is indeed a simplex, then we should have 
                # found exactly s edges 
                if SubedgeCount == s:
                    self.simplexCounts[s] = self.simplexCounts[s] + 1
            
        # Every simplex is counted 4 times so we need to divide by 4
        self.simplexCounts[s] = self.simplexCounts[s] // 4

    def simplexApproxBadVariance(self, s, e, d):
        """
        (e,d) approximation to the number of s-simplices in the hypergraph using the 
        method with the poor variance

        Args
            s: size of simplex being counted
            e: precision
            d: confidence
        """
        # By the median of means theorem, we need to take the means of running our 
        # estimator algorithm for a number of repetitions and then take the median 
        # of those means over a set of repetitions
        # print(self.m)
        repetitionsMeans = math.ceil(self.m**(2-1/(s-1))/(e*e*self.simplexCounts[s]))
        repetitionsMedian = math.ceil(math.log(1/d))
        # print(repetitionsMeans, repetitionsMedian)
        meansList = []
        for i in range(repetitionsMedian):
            meansSingle = 0
            print("Median of means: ", i,"-th means calculation")
            for j in range(repetitionsMeans):
                meansSingle = meansSingle + self.simplexApproxBadVarianceSingle(s, False)
            meansSingle = meansSingle/repetitionsMeans
            meansList.append(meansSingle)

        # How much memory does this method use?
        self.badVarianceMemoryUsage = repetitionsMeans*repetitionsMedian
        print("The method with bad variance uses approximately ", self.badVarianceMemoryUsage, "words of memory.")
        
        return statistics.median(meansList)


    def reservoirSampling (self):
        """
        Basic reservoir sampling method in streaming
        
        Returns:
            a uniformly sampled item from the stream of hyperedges
            the number of elements in the stream
        """
        # Basic stream index running through the stream
        idx = 1
        # The edge we will sample
        chosenEdge = []

        # Reservoir sampling requires only one pass through the input stream
        self.file_reader.seek(0, 0)
        edgeString = self.file_reader.readline()
        while edgeString:
            edge = edgeString.split()
            # Retain edge with probability 1/idx
            coinToss = random.random()
            if coinToss < 1/idx:
                chosenEdge = copy.deepcopy(edge)
            idx = idx+1
            edgeString = self.file_reader.readline()

        return (chosenEdge, idx-1)

    def simplexApproxBadVarianceSingle (self, s, showProgress):
        """
        Approximately Count the simplices of the streamed hypergraph,
        using the Bad Variance estimator algorithm

        Args
            s: the simplex size to be counted
            showProgress: should we show the progress?
        """

        if showProgress:
            print("BAD VARIANCE METHOD: Approximating simplex count of size ",s)
        # s is outside of range of simplex we count
        if s <= 2 or s >= self.mhs + 3:
            return
        
        # Pass 1: pick an edge using reservoir sampling
        if showProgress:
            print("Pass 1: pick an edge using reservoir sampling")
        (e,m) = self.reservoirSampling()
        if showProgress:
            print("Number of edges = ",m)

        # Pass 2: Calculate all the degrees in the chosen edge
        if showProgress:
            print("Pass 2: Calculate all the degrees in the chosen edge")
        degreesInEdge = [0 for i in range(len(e))]
        self.file_reader.seek(0, 0)
        edgeString = self.file_reader.readline()
        while edgeString:
            edge = edgeString.split()
            for i in range(len(e)):
                if e[i] in edge:
                    degreesInEdge[i] = degreesInEdge[i] + 1
            edgeString = self.file_reader.readline()

        # Pass 3: 
        if showProgress:
            print("Pass 3")
        # r is the variance-controlling parameter
        r = math.ceil(min(degreesInEdge)/math.sqrt(m))
        if showProgress:
            print ("r = ", r)
        # Re-arranging the vertices in e so that they are in increasing order of degree
        e = [x for _,x in sorted(zip(degreesInEdge,e))]
        degreesInEdge.sort()
        
        # e[0]-Neighbor samples to be populated
        X = [0 for k in range(r)]
        # The idx we will use for the reservoir sampling
        neighborhoodIdx = 1
        # Resetting the file reader to the beginning of the file
        self.file_reader.seek(0, 0)

        
        progress = 0
        percentProgress = 0
        edgeString = self.file_reader.readline()
        while edgeString:
             # Keep track of progress
            progressed = False
            while progress >= percentProgress*m:
                percentProgress = percentProgress + 0.01
                progressed = True
            if progressed and (int(percentProgress*100))%20 == 0 and showProgress:
                print("Progress: ", int(percentProgress*100), "%")
            progress = progress + 1

            edge = edgeString.split()
            edgeString = self.file_reader.readline()
            # only look at neighbors of e[0]
            if e[0] not in edge:
                continue
            # form a stream of the neighbors of e[0] that are not in e
            edge.remove(e[0])
            edge = [v for v in edge if not v in e]
            # sample a neighbor from this virtual stream uniformly at random
            # using reservoir sampling
            for v in edge:
                for k in range(r):
                    coinToss = random.random()
                    if coinToss < 1/neighborhoodIdx:
                        X[k] = v
                neighborhoodIdx = neighborhoodIdx+1
        
        # Save the virtual stream length (number of neighbors of e[0] that are not in e)
        # This will be used for the estimator
        neighborhoodStreamLength = neighborhoodIdx - 1

        # Pass 4: Check if the assembled simplex is indeed a valid simplex
        if showProgress:
            print("Pass 4: Check if the assembled simplex is indeed a valid simplex")

        # As in the naive algorithm, we will be counting how many edges are present in a simplex
        simplexPresentEdges = [0 for k in range(r)]
        # Keep track of the simplices we will be trying to detect
        proposedSimplices = []
        for k in range(r):
            proposedSimplex = copy.deepcopy(e)
            proposedSimplex.append(X[k])
            proposedSimplices.append(proposedSimplex)
        if showProgress:
            print("Proposed Simplices: ", proposedSimplices)
        # We will also need to keep track of the degrees of the X-vertices
        degreesOfXs = [0 for k in range(r)]
        # ... and also keep track of their e[0] codegrees
        codegreesOfXs = [0 for k in range(r)]

        progress = 0
        percentProgress = 0
        self.file_reader.seek(0, 0)
        edgeString = self.file_reader.readline()
        while edgeString:
            # Keep track of progress
            progressed = False
            while progress >= percentProgress*m:
                percentProgress = percentProgress + 0.01
                progressed = True
            if progressed and (int(percentProgress*100))%20 == 0 and showProgress:
                print("Progress: ", int(percentProgress*100), "%")
            progress = progress + 1

            edge = edgeString.split()
            edgeString = self.file_reader.readline()
            # For each proposed simplex, check if the incoming edge is a part of it
            for k in range(r):
                if set(edge).issubset(set(proposedSimplices[k])):
                    simplexPresentEdges[k] = simplexPresentEdges[k] + 1
                # We also need to calculate the degrees of the X[k]s
                if X[k] in edge:
                    degreesOfXs[k] = degreesOfXs[k] + 1
                    # We also need to calculate the e[0]-codegrees of the X[k]s
                    if e[0] in edge:
                        codegreesOfXs[k] = codegreesOfXs[k] + 1
        
        # After watching the whole stream pass by, we need to check the following conditions
        # for each k = 1,2,...,r:
        #   1. X[k] > e[-1] in the total ordering based on degree
        #   2. The k-th simplex is complete
        # Based on those conditions, we set the value of our estimators
        Z = [0 for k in range(r)]
        for k in range(r):
            if ((degreesOfXs[k] > degreesInEdge[-1]) or (degreesOfXs[k]==degreesInEdge[-1] and X[k]>e[-1])) and simplexPresentEdges[k] == s:
                if showProgress:
                    print("Simplex detected")
                Z[k] = neighborhoodStreamLength / codegreesOfXs[k]
        
        # Take the average - this is our final estimator
        if showProgress:
            print("---------------------------------------")
        return (m*sum(Z))/len(Z)

try:
    generator = HypergraphGenerator(30, 0.7)
    m = generator.generate("C:/Users/themi/OneDrive/Desktop/Thesis/Code/testing_hypergraph_generation_1.txt")
    stream = Stream("C:/Users/themi/OneDrive/Desktop/Thesis/Code/testing_hypergraph_generation_1.txt",
                    30,
                    m,
                    3)
    stream.CountAllSimplices(4)
    print("Number of 4-simplices = ", stream.simplexCounts[4])
    print("Approximation with bad variance = ", stream.simplexApproxBadVariance(4, 0.1, 0.01))

finally:
    stream.done()
