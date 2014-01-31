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
        return ''.join(`self.seq`)[::(-2)*(-self.orientation)-1]

    def reverse(self):
        self.orientation = not self.orientation


class Chromosome:
    def __init__(self, nbGen, sizeGen, sdGen, alphabet='ATCG'): #%ATGC
        self.alphabet=np.array(list(alphabet))
        self.nbGen=nbGen
        self.sizeGen=sizeGen
        self.sdGen=sdGen
        self.chromGeneration()
        self.thisChromSize=np.sum(len(gen.seq) for gen in self.seq)
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
            
    def mutation(self,pMut):
        i=0
        nbMut=self.thisChromSize*pMut
        while i<nbMut: 
            gIndex=rand.randint(0,self.nbGen) #index of the gene to mutate
            mutGen=self.seq[gIndex].seq #sequence of the gene to mutate
            
            nIndex=rand.randint(0,len(mutGen)) #indice of the nucleotide to mutate

            mutGen[nIndex] = (mutGen[nIndex]+rand.randint(1,4))%4 #mutate randomly to one of the other bases  		
            i+=1
        return int(nbMut+1) #actual number of mutations processed


if __name__=='__main__':
    
    # Simulate the indicated number of inversion, on a genome containing the indicated number of genes, and measure the estimated distance between the genom and his ancestor every "measureStep" inversions.
    #Indication : 100 measures on a 10000 genes genomes last for around 25 minutes
    #Result is stored in a file with the specified fileName, one measure by row.

    def simulEvolInversions(nbGenes, nbInversions, measureStep, fileName):
        def rearrange(a,b,chrom,graph):
            graph.rearrange(chrom.seq[a],chrom.seq[b-1],chrom)
            chrom.rearrange(a,b)

        nbGen=nbGenes
        chrom=Chromosome(nbGen,5,1)
    #    print chrom
        gr=Graph(chrom)
    #    print gr
        rangeNbGen=range(nbGen)
        f=open(fileName,"w")
        for i in range(nbInversions):
            a,b=random.sample(rangeNbGen,2)
            rearrange(min(a,b),max(a,b),chrom,gr)
            if i%measureStep==0:
                #        print i
                f.write("%d\n"%gr.distance())
        f.close()

    # Simulate the indicated number of generations, processing random mutations (with the mutation rate pMut) on a genome containing the indicated number of genes,
#with the indicated number of bases in each gene (taken in a normal low of sd "sdGen")
# and measure the estimated distance between the genom and his ancestor every "measureStep" generations.
    #Result is stored in a file with the specified fileName, on each row you'll find : nbGenerations - realNumber of mutations - observed number of mutations.

    def simulEvolMutations(nbGen,sizeGen,sdGen,pMut,nbGenerations,measureStep,fileName):
        def alignement(chrom1,chrom2):
            align=[gen1.seq-gen2.seq for gen1,gen2 in zip(chrom1.seq,chrom2.seq)]
            return np.sum([len(np.flatnonzero(seq)) for seq in align] )

        originalChrom=Chromosome(nbGen,sizeGen,sdGen,alphabet=(0,1,2,3)) # 1='A' , 2='C', 3='G', 0='T'
        chromMut=copy.deepcopy(originalChrom)

        f=open(fileName,"w")
        nReel=0
        for i in range(nbGenerations):
            nReel+=chromMut.mutation(pMut)
            if i%measureStep==0:
                f.write(str(i)+" "+str(nReel)+"  "+str(alignement(originalChrom,chromMut))+"\n")
        f.close()
    
    
    #Simulations
    simulEvolInversions(nbGenes=10,nbInversions=10,measureStep=1,fileName="dist.data")
    simulEvolMutations(nbGen=100,sizeGen=100,sdGen=10,pMut=1e-4,nbGenerations=1000,measureStep=10,fileName='mut.data')
