import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import cobra.io
from numpy import *
from os.path import join
from scipy.stats import *
from bisect import bisect_left




######Essential Gene Discovery Functions#########
def binary_search(a, x, lo=0, hi=None):  # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)
    pos = bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1)  # don't walk off the end  #https://stackoverflow.com/questions/212358/binary-search-bisection-in-python
"""
Creates Essential Gene Data for CS models. Saves the results as a text file named EssGeneACCESSIONNUMBER_SAMPLENUMBER.txt with the 
deleted gene name in the first column and the flux through biomass on the second column

Inputs: assent - GEO accession number that corressponds to the name on the .mat file
modelType - (optional) type of model file to be read (now only supports .mat)
numSamples - # of samples that correspond to the accession number
path - File path for location of the .mat models
solver - optimizer to be used for doing the simulation 'gurobi is the default'

Outputs: None
"""
def createEssentialGeneDataAssent(assent,modelType = '.mat',numSamples =1 ,path=None,solver='gurobi'):
    newFile = open(('EssGene'+assent+'.txt'),'w')
    for x in [x+1 for x in range(numSamples)]:
        fileName = join(path,assent+ '_' + str(x) + modelType)
        CSmodel = cobra.io.load_matlab_model(fileName)
        CSmodel.solver = solver   #load model
        res = cobra.flux_analysis.single_gene_deletion(CSmodel) #perform single gene deletion simulation
        res.sort_index()  #sort the results of the simulation
        results = []
        for y in res.to_records(index=True):
            results.append(tuple(y)) #gather geneName,flux tuples
        newFile.write('\n-----------------------\n'+ assent+'_'+str(x)+'\n')#write sample label
        [newFile.write(str(x[0]) + ' ' + str(x[1])+'\n') for x in results] #write results to file

    newFile.write('\n-----------------------\n')#write final marker
    newFile.close()

"""
No longer functional. Original Purpose is to find the number of changed genes given a set of unique genes. Not compatable with flux variatios=ns
def findNumberofHits(uniqueData):
    pings = []
    for x in uniqueData:
        for r in x:
            for g in x[r]:
                count = 0
                if len(pings) == 0:
                    pings = [GeneHits(g)]
                else:
                    for j in pings:
                        if j.name == g:
                            count = count + 1
                    if count == 0:
                        pings.append(GeneHits(g))
                    else:
                        for j in pings:
                            if j.name == g:
                                j.addToCount()
    return pings

"""


"""
Function that calculates the results of a single deletion simulation for a cobra model. Creates a text file in the same format as above.
Inputs: model - cobra model for which the simulation should be performed
label - label for the resulting file output i.e. EssGeneLABEL.txt
solver - solver for use in simulation default is gurobi
"""
def createEssentialGeneModel(model,label,solver = 'gurobi'):
    model.solver = solver
    res = cobra.flux_analysis.single_gene_deletion(model)#perform simulation
    newFile = open(('EssGene' + label + '.txt'), 'w')
    res.sort_index()
    results = []
    for x in res.to_records(index=True):
        results.append(tuple(x)) #gather geneName flux level tuple
    [newFile.write('\n' + str(x[0]) + ' ' + str(x[1])) for x in results] #write results
    newFile.write('\n')
    newFile.close()

"""   
def createEssentialGeneModelDD(model,label,type = 'smbl',solver = 'gurobi'):
    model.solver = solver
    res = cobra.flux_analysis.double_gene_deletion(model)
    newFile = open(('EssGene'+label+'.txt'),'w')
    res.sort_index()fdsafd
    results = []
"""

    ##### Analysis Tools #########
"""   
Class the handles the information from simulation results. Organizes flux levels for a series of hits. Allows for consolidated information and output
"""
class GeneHits:
    def __init__(self , differences,  name):
        self.name = name
        self.count = 1
        self.differences = [differences]
    def addToCount(self):
        self.count = self.count + 1
    """ Returns a formatted string containing the name of the gene, number of hits, mean, std, sum, and pvalue using one sample t test"""
    def letsPrint(self):
        if len(self.differences) > 1 and not(mean(self.differences) == 0.0):    #TODO add in null data to the differences array to show more accurate statistic; multiply by control coeff pad with 0
            _,pval = ttest_1samp(self.differences,0.0)
        else:
            pval = 1.0
        return self.name + ': %d' % self.count + '; mean: %.12lf ; std: %.12lf ; sum: %.12lf ; pvalue %.08f' %(mean(self.differences),std(self.differences),sum([abs(x) for x in self.differences]),pval)
    def areEqual(self,x):
        if x.name == self.name:
            return 1
        else:
            return 0
    """Combines two GeneHits objects"""
    def combine(self,x):
        self.addToCount()
        self.differences = self.differences + x.differences

"""
Class for the characterization of essential genes for a general metabolic model. Processes the data in the EssGene file. 
Data is acessible through the getData() function wher the genes and corresponding fluxes are stored as a dictionary
"""
class ComparisionGene:
    def __init__(self,label,information):
        self.normalFile = open('EssGene'+label+'.txt', 'r')
        self.normalData = self.normalFile.readlines()
        try:
            while(True):
                self.normalData.remove('\n') #clean up file
        except:
            pass
        for x in range(len(self.normalData)):
                self.normalData[x] = self.normalData[x][:-1] #remove final characters
        self.normalData = [x.split(" ") for x in self.normalData]#separate name and values
        self.normalData = {x[0]:float(x[1]) for x in self.normalData}#creat dictionary
        self.normalFile.close()
        self.information = information #stores other string information that corresponds to the general model.
        self.normalFile.close()
    def getData(self):
        return self.normalData
"""
#Class for comparision and analysis of essential gene data from a context specific model. Processes the flux levels in the text file
. Allows for comparision of data from the file with a general model data. Can handle the control and experimental data

"""
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
                    tempData2 = [tempData[0]]
                    for y in [x+1 for x in range(len(tempData)-1)]:
                        spliter = tempData[y].split(" ")
                        tempData2.append([spliter[0],float(spliter[1])])
                    self.sampleEG.append(tempData2)
                tempData = list()
            else:
                tempData.append(x)
        try:
            self.sampleEG.remove([])
        except:
            pass
        self.sampleEG = {x[0]:x[1:] for x in self.sampleEG}
        for x in self.sampleEG:
            self.sampleEG[x] = {y[0]:y[1] for y in self.sampleEG[x]}
        keysC = []
        keysE = []

        for i in range(len(self.control)):
            if self.control[i] == 0:
                keysC.append(self.accession + '_' + str(i+1))
            else:
                keysE.append(self.accession+'_'+str(i+1))
        #totList =[]
        #[totList.append(set(self.sampleEG[x])) for x in keysC]
        #self.controlData = [list(x) for x in totList[:]]
        #totList = list()
        #[totList.append(set(self.sampleEG[x])) for x in keysE]
        #self.experimentalData = [list(x) for x in totList[:]]
        self.controlData = [self.sampleEG[x] for x in keysC] #TODO need to handle all control or all experimental samples (index out of bounds)
        self.experimentalData = [self.sampleEG[x] for x in keysE]


    def findUniqueFrom(self,geneTot):
        uniqueGene = []
        controlUnique = []
        experimentalUnique = []
        i = 0
        for x in self.sampleEG:
            unique = list(set(self.sampleEG[x])-set(geneTot))
            uniqueGene.append([x,unique[:]])
            if self.control[i] == 1:
                controlUnique.append([x,unique[:]])
            else:
                experimentalUnique.append([x,unique[:]])
            i+=1
        totalUnique = {x[0]:x[1:][0] for x in uniqueGene}
        controlUnique = {x[0]:x[1:][0] for x in controlUnique}
        experimentalUnique = {x[0]:x[1:][0] for x in experimentalUnique}
        return totalUnique,controlUnique,experimentalUnique

    def findChangedFluxGenes(self,generalGene,group):
        results = []
        totalResults = []
        data = self.sampleEG
        i = 0
        for x in data:
            coeff = self.control[i]
            for y in data[x]:
                if not(data[x][y] == generalGene[y]):
                    if len(results) == 0:
                        results.append(GeneHits(data[x][y] - generalGene[y],y))
                    elif coeff == group:
                        count = 0
                        for z in results:
                            if z.name == y:
                                z.combine(GeneHits((data[x][y] - generalGene[y]),y))
                                count +=1
                        if count == 0:
                            results.append(GeneHits((data[x][y] - generalGene[y]),y))
            i +=1
        numE = 0
        for x in self.control:
            if x == group:
                numE += 1
        for x in results:
            if len(x.differences) > numE:
                temp = array(x.differences)
                temp.sort()
                temp = list(temp)
                x.differences = temp[0:numE-1]
            while len(x.differences) < numE:
                x.differences.append(0.0)

        return results,numE


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












