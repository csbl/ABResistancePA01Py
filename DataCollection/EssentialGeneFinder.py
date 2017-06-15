import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import cobra.io
from numpy import *
import os
from os.path import join


from bisect import bisect_left
######Essential Gene Discovery Tools#########
def binary_search(a, x, lo=0, hi=None):  # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)
    pos = bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1)  # don't walk off the end  #https://stackoverflow.com/questions/212358/binary-search-bisection-in-python

def createEssentialGeneDataAssent(assent,modelType = '.mat',numSamples =1 ,path=None,solver='gurobi'):
    EssGeneDataCollection = []
    for x in [x+1 for x in range(numSamples)]:
        fileName = join(path,assent+ '_' + str(x) + modelType)

        CSmodel = cobra.io.load_matlab_model(fileName)

        CSmodel.solver = solver

        res = cobra.flux_analysis.single_gene_deletion(CSmodel)

        essential = res['flux'] <= 0.001 #Handle numerical error

        essGenes = res[essential]

        essGenes = sort(array(essGenes.index.tolist()))

        EssGeneDataCollection = EssGeneDataCollection + ['\n' + assent + str(x) + '\n-----------------------\n'] +  essGenes.tolist()

    newFile = open(('EssGene'+assent[:-1]+'.txt'),'w')

    [newFile.write('\n'+x) for x in EssGeneDataCollection]

    newFile.close()

def createEssentialGeneModel(model,label,type = 'smbl',solver = 'gurobi'):
    model.solver = solver
    res = cobra.flux_analysis.single_gene_deletion(model)
    essential = res['flux'] == 0
    essGenes = res[essential]
    essGenes = sort(array(essGenes.index.tolist()))
    newFile = open(('EssGene' + label + '.txt'), 'w')
    [newFile.write('\n' + x) for x in essGenes]
    newFile.close()
    ##### Analysis Tools #########
class GeneHits:
    def __init__(self , name):
        self.name = name
        self.count = 1
       # self.index = name[:6]
    def addToCount(self):
        self.count = self.count + 1
    def letsPrint(self):
        print self.name + ': %d' % self.count
class ComparisionGene:
    def __init__(self,label,information):
        self.normalFile = open('EssGene'+label+'.txt', 'r')
        self.normalData = self.normalFile.readlines()
        for x in self.normalData:
            if x == '\n' or x.startswith('-') == 1:
                self.normalData.remove(x)
        self.normalFile.close()
        self.information = information
    def getData(self):
        return self.normalData
class CSComparisionGene:
    def __init__(self,assention,controlIndex,information):
        self.assention = assention
        self.information = information
        self.file = open('EssGene'+assention+'.txt','r')
        self.dataTot = self.file.readlines()
        self.sampleData = []
        self.controlIndex = controlIndex
        tempData = []
        for x in self.dataTot:
            if x == '\n':
                self.dataTot.remove(x)
            elif x.startswith('G') == 1:
                tempData.append(x)
            elif x.startswith('-') == 1:
                self.sampleData.append(tempData)
                tempData = []
            else:
                tempData.append(x[:-1])
        self.sampleData = {x[0]:x[1:] for x in self.sampleData}
