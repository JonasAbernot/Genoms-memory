#-*- coding:utf-8 -*-
import numpy as np
import copy

class Graph:
    def __init__(self,chrom):
        self.adjacencyMatrix = self.chromToGraph(chrom)
        self.vertRange=range(len(self.adjacencyMatrix))

    def __repr__(self):
        return str(self.adjacencyMatrix)

    #Build the INITIAL arrangement graph, with all genes oriented the same way
    #Matrix of the shape [0, 0, 0, ..., 0, 0, 0]
    #                    [0, 0, 2, ..., 0, 0, 0]
    #                    [0, 2, 0, ..., 0, 0, 0]
    #                     .                   .
    #                     .                   . etc
    def chromToGraph(self,chrom): 
        n=chrom.nbGen
        adjMat=np.zeros((2*n,2*n),dtype=np.int8)
        for i in range(1,2*n-1,2):
            adjMat[i,i+1]=2
            adjMat[i+1,i]=2
        return adjMat

    #Adapt the adjacency matrix when the genome undergoes a inversion
    def rearrange(self,a,b,chrom):

        #find the vertex wich share an edge with the specified one in the nowdays arrangement, if it exists
        def findLinked(a,adjVec):
            alinks=np.flatnonzero(adjVec).tolist() # all neighbourgs of vertex a
            a2=None
            if len(alinks)>1:# if vertex a is linked to several (i.e. two) vertices
                alinks.remove(a+2*(a%2)-1) # Don't take account of the one from the ancestral genom
                a2=alinks[0]
            elif len(alinks)==1:# and if it's linked to only one, 
                a2=alinks[0]
                if adjVec[a2]==1 and a2==(a+2*(a%2)-1):#check it's not the link from the ancestral genom
                    a2=None
            return a2

        #step 1: find vertices concerned by the inversion
        a1 = 2*a.id + (not a.orientation)
        b1 = 2*b.id + b.orientation

        #step 2:find vertices linked to the first ones in the nowdays arrangement
        a2=findLinked(a1,self.adjacencyMatrix[a1])
        b2=findLinked(b1,self.adjacencyMatrix[b1])

        #step 3:cross the edges
        if a2 != None:
            self.adjacencyMatrix[a1,a2]-=1
            self.adjacencyMatrix[a2,a1]-=1
            self.adjacencyMatrix[b1,a2]+=1
            self.adjacencyMatrix[a2,b1]+=1
        if b2 != None:
            self.adjacencyMatrix[b1,b2]-=1
            self.adjacencyMatrix[b2,b1]-=1
            self.adjacencyMatrix[a1,b2]+=1
            self.adjacencyMatrix[b2,a1]+=1

            
    #return a list of connected sub-graphs of the graph
    def groupAdjacencyMatrix(self):

        #return the maximal connected subgraphs the vertex i belongs to (BFS)
        def buildGroup(self,i):
            group=set([i])
            toAddToGroup=set(np.flatnonzero(self.adjacencyMatrix[i]))

            while any([j not in group for j in toAddToGroup]):
                toAddToGroup = toAddToGroup - (toAddToGroup & group) 
                group = group | toAddToGroup
                newToAdd = set([])

                for j in toAddToGroup:
                    newToAdd = newToAdd | set([k for k in np.flatnonzero(self.adjacencyMatrix[j])])

                toAddToGroup = newToAdd

            return group
        
        connectedPaths=[]
        explored=[False]*len(self.adjacencyMatrix)

        for i in self.vertRange:
            if not explored[i]:
                group=buildGroup(self,i)
                for j in group:
                    explored[j]=True
                connectedPaths.append(group)

        return connectedPaths

    #compute the minimum number of inversions necessary to come back to the ancestral shape of the genome
    def distance(self):
        #d=n-(c+i/2)

        #initialisation:
        n=len(self.adjacencyMatrix)/2 #number of genes
        c=0 #number of cycles
        i=0 #nombre of uneven paths

        groups = self.groupAdjacencyMatrix()

        #thanks to the special structure of the graph, we know the subgraph can only be either a simple unbranched tree, including a single node, or a cycle.
        for group in groups:
            if len(group) == 1:
                i+=1
            else:
                sums=[np.sum(self.adjacencyMatrix[j]) for j in group] #number of degrees of each vertex from the connected subgraph
                if 1 in sums: #one of the vertex is of degree 1, so the subgraph is not a cycle
                    if len(sums)%2==1:
                        i+=1
                else: #it's a cycle.
                    c+=1

        return n-(c+i/2)
        
