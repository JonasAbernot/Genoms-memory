#-*- coding:utf-8 -*-

import numpy.random as rand
import numpy as np
import copy
import random
#___________
from rearrangementGraph import Graph

class Gene:
    def __init__(self,seq,id):
        self.seq=seq
        self.orientation=True
        self.id=id

    def __repr__(self):
        return ''.join(self.seq)[::(-2)*(-self.orientation)-1]

    def reverse(self):
        self.orientation = not self.orientation


class Chromosome:
    def __init__(self, nbGen, sizeGen, sdGen, alphabet='ATCG'): #%ATGC
        self.alphabet=np.array(list(alphabet))
        self.nbGen=nbGen
        self.sizeGen=sizeGen
        self.sdGen=sdGen
        #also contain the genes list (self.seq)

    def genGeneration(self,id):
        thisGenSize=rand.normal(scale = self.sdGen,
                                loc = self.sizeGen, 
                                size = 1) 
        return Gene(self.alphabet[rand.randint(4,size=thisGenSize)],id)
		
		    
    def chromGeneration(self):
        self.seq=[self.genGeneration(i) for i in range(self.nbGen)]

    def __repr__(self):
        return `self.seq`

    def rearrange(self,a,b):
        for gene in self.seq[a:b]:
            gene.reverse()
        if a!=0:
            self.seq = self.seq[:a] + self.seq[b-1:a-1:-1] + self.seq[b:]
        else: 
            self.seq = self.seq[:a] + self.seq[b-1::-1] + self.seq[b:]
            
if __name__=='__main__':
    def rearrange(a,b,chrom,graph):
        graph.rearrange(chrom.seq[a],chrom.seq[b-1],chrom)
        chrom.rearrange(a,b)

    nbGen=10
    chrom=Chromosome(nbGen,5,1)
    chrom.chromGeneration()
#    print chrom
    gr=Graph(chrom)
#    print gr
    rangeNbGen=range(nbGen)
    f=open("dist.data","w")
    for i in range(10*nbGen):
        a,b=random.sample(rangeNbGen,2)
        rearrange(min(a,b),max(a,b),chrom,gr)
        if i%1==0:
    #        print i
            f.write("%d\n"%gr.distance())
    f.close()
