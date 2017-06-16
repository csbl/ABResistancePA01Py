import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import cobra.io
from numpy import *
from os.path import join


from bisect import bisect_left
######Essential Gene Discovery Tools#########
def binary_search(a, x, lo=0, hi=None):  # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)
    pos = bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1)  # don't walk off the end  #https://stackoverflow.com/questions/212358/binary-search-bisection-in-python

def createEssentialGeneDataAssent(assent,modelType = '.mat',numSamples =1 ,path=None,solver='gurobi'):
    essGeneDataCollection = []
    for x in [x+1 for x in range(numSamples)]:
        fileName = join(path,assent+ '_' + str(x) + modelType)
        CSmodel = cobra.io.load_matlab_model(fileName)
        CSmodel.solver = solver
        res = cobra.flux_analysis.single_gene_deletion(CSmodel)
        essential = res['flux'] <= 0.001 #Handle numerical error
        essGenes = res[essential]
        essGenes = sort(array(essGenes.index.tolist()))
        essGeneDataCollection = essGeneDataCollection + ['\n-----------------------\n'+ assent+'_'+str(x)] +  essGenes.tolist()
    essGeneDataCollection = essGeneDataCollection + ['\n-----------------------\n']

    newFile = open(('EssGene'+assent+'.txt'),'w')

    [newFile.write('\n'+x) for x in essGeneDataCollection]

    newFile.close()

def createEssentialGeneModel(model,label,type = 'smbl',solver = 'gurobi'):
    model.solver = solver
    res = cobra.flux_analysis.single_gene_deletion(model)
    essential = res['flux'] <= 0.001
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

#
#Class for the characterization of essential genes for a general metabolic model. Processes the data in the EssGene file.
#
class ComparisionGene:
    def __init__(self,label,information):
        self.normalFile = open('EssGene'+label+'.txt', 'r')
        self.normalData = self.normalFile.readlines()
        try:
            while(True):
                self.normalData.remove('\n')
        except:
            pass
        for x in range(len(self.normalData)):
                self.normalData[x] = self.normalData[x][:-1]
        self.normalFile.close()
        self.information = information
        self.normalFile.close()
    def getData(self):
        return self.normalData

#Class for comparision and analysis of essential gene data from a context specific model
#
#
#

class CSGene:
    def __init__(self,accession,control,media = 'LB',information='NO INFORMATION'): #TODO remove default value for media
        self.accession = accession
        self.information = information
        self.file = open('EssGene'+self.accession+'.txt','r')
        self.dataTot = self.file.readlines()
        self.control = control
        self.media = media
        try:
            while(True):
                self.dataTot.remove('\n')
        except:
            pass
        for x in range(len(self.dataTot)):
            self.dataTot[x] = self.dataTot[x][:-1]
        self.file.close()
    def processSamples(self):
        self.sampleEG = []
        tempData = []
        for x in self.dataTot:
            if x[0] == '-':
                if len(self.sampleEG) == 0:
                    self.sampleEG = [tempData]
                else:
                    self.sampleEG.append(tempData)
                tempData = list()
            else:
                tempData.append(x)
        self.sampleEG.remove([])
        self.sampleEG = {x[0]:x[1:] for x in self.sampleEG}
        keysC = []
        keysE = []
        for i in range(len(self.control)):
            if self.control[i] == 1:
                keysC.append(self.accession + '_' + str(i+1))
            else:
                keysE.append(self.accession+'_'+str(i+1))
        totList = []
        [totList.append(set(self.sampleEG[x])) for x in keysC]
        self.controlData = [list(x) for x in totList[:]]
        totList = list()
        [totList.append(set(self.sampleEG[x])) for x in keysE]
        self.experimentalData = [list(x) for x in totList[:]]


    def findUniqueFrom(self,geneTot):
        uniqueGene = []
        controlUnique = []
        experimentalUnique = []
        i = 0
        for x in self.sampleEG:
            unique = list(set(geneTot) -  set(self.sampleEG[x]))
            uniqueGene.append([x,unique])
            if self.control[i] == 1:
                controlUnique.append([x,unique])
            else:
                experimentalUnique.append([x,unique])
            i+=1
        totalUnique = {x[0]:x[1:][0] for x in uniqueGene}
        controlUnique = {x[0]:x[1:][0] for x in controlUnique}
        experimentalUnique = {x[0]:x[1:][0] for x in experimentalUnique}
        return totalUnique,controlUnique,experimentalUnique

    def findUniqueFromControl(self):
        try:
            totList = []
            [totList.append(set(x)) for x in self.controlData]
            totControl = list(set.intersection(*totList)) #find intersection of all genes for the control
            totList = list()
            [totList.append(set(x)) for x in self.experimentalData]
            interExp = list(set.intersection(*totList)) #find genes that are common to all experimental sample essential genes
            return list(set(interExp) - set(totControl)) #find members of the shared genes among the intersection that are not in the control genes
        except :
            print 'Error: You must run CSGene.processSamples() first'
            return -1


#def getExperimentalSamples():










