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
            alinks=np.flatnonzero(adjVec).tolist()
            a2=None
            if len(alinks)>1:
                alinks.remove(a+2*(a%2)-1)
                a2=alinks[0]
            elif len(alinks)==1:
                a2=alinks[0]
                if adjVec[a2]==1 and a2==(a+2*(a%2)-1):
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

        #return the maximal connected subgraphs the vertex i belongs to
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

        for group in groups:
            if len(group) == 1:
                i+=1
            else:
                sums=[np.sum(self.adjacencyMatrix[j]) for j in group]
                if 1 in sums:
                    if len(sums)%2==1:
                        i+=1
                else:
                    c+=1

        return n-(c+i/2)
        
