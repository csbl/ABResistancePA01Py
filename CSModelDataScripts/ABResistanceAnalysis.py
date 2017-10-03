import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import cobra.io
from numpy import *
from os.path import join
from scipy.stats import *
from bisect import bisect_left
import cobra.manipulation
from copy import deepcopy
import scipy.stats

######Essential Utility Functions#########
def binary_search(a, x, lo=0, hi=None):  # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)
    pos = bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1)  # don't walk off the end  #https://stackoverflow.com/questions/212358/binary-search-bisection-in-python
"""Find the 2 side pvalue for a 2 sample ttest with unequal variances given the means (y1,y2), number of samples (n1,n2), and variances (s1,s2)"""
def t_testUnpaired_fromSum(y1,s1,n1,y2,s2=-1,n2= -1,):
    if s2 == -1: s2 = s1
    if n2 == -1: n2 = n1
    T = (y1-y2)/sqrt(s1/n1 + s2 / n2)
    v = ((s1/n1+s2/n2)**2)/(((s1/n1)**2)/(n1-1)+((s2/n2)**2)/(n2-1))
    return scipy.stats.t.sf(abs(T),v)[0]*2
#####Simulation and Analysis Function for Context Specific Models###########
"""Creates consensus model by deleting gene list prescribed in geneList, or in typeGenes.txt which must be created. 
 The model is optimized using FBA, the removed reactions are tabulated, single gene and reaction deletion simulations are performed,
as well as FVA. The results are all saved in text or csv files which are analyzed in the ModelComparison.py file. """
def consensusModelCreator(model,type,geneList = [1]):
    if geneList[0] == 1:
        geneListF = open(type+'Genes.txt','r')
        geneList = geneListF.readlines()
        geneList = [x[:-1] for x in geneList]
        geneListF.close()
        print "Making standard Model\n"
    #print geneList
    cobra.manipulation.delete_model_genes(model,geneList)#create consensus model
    model.optimize()# peform FBA

    removedReactionsID = []
    for x in model.reactions:
        if x.bounds == (0,0):
            removedReactionsID.append(x.id) #tabulate removed reactions
    outputFile = open(type +'ReactionsKO.txt','w')
    outputFile.write('ReactionID\n')
    print removedReactionsID
    [outputFile.write(x+'\n') for x in removedReactionsID]
    outputFile.close()

    res = cobra.flux_analysis.single_gene_deletion(model)#perform single gene deletion
    res.to_csv(type+'GeneDeletion.csv')
    res = cobra.flux_analysis.single_reaction_deletion(model)#perform single reaction deletion
    res.to_csv(type + 'RxnDeletion.csv')
    print "done"
    res = cobra.flux_analysis.flux_variability_analysis(model,loopless = True)#perform FVA
    res.to_csv(type + 'FVA.csv')
"""Finds deleted reactions that are shared across all media conditions for the control and resistant groups
"""
def creationOfDeletedReaction():
    assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
    reference = {'GSE90620': [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1], 'GDS3572': [1, 1, 1, 0, 0, 0],
                 'GSE65870': [1, 1, 1, 1, 1, 1], 'GDS4244': [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                 'GSE30021': [0, 0, 0, 0, 0, 0, 1, 1, 1]}
    referenceMed = {'GSE90620': 1, 'GDS3572': 0, 'GSE65870': 2, 'GDS4244': 1, 'GSE30021': 1}#media code: 1=LB, 2=PBM,0 = MH
    nsamples = [12,12,6,9,6]
    deletions = [[set(),set()],[set(),set(),set()]]
    model = cobra.io.read_sbml_model("iPAE1146.xml")
    for i in range(len(assent)):#iterate through each sample
        res = findDeletedReactions(assent[i],model,nsamples[i],reference[assent[i]],'gurobi','C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')#find deleted reactions for the sample set
        if not len(res[0]) == 0 : deletions[0][referenceMed[assent[i]]] = deletions[0][referenceMed[assent[i]]].union(res[0])#add on unique values to control
        else: print "None found in control: %d" % i
        if not len(res[1]) == 0 : deletions[1][referenceMed[assent[i]]] = deletions[1][referenceMed[assent[i]]].union(res[1])#add on unique values to experimental
        else: print "None found in exp: %d" % i

    return [deletions[0][1] & deletions[0][1],deletions[1][0] & deletions[1][1] & deletions[1][2]] #find shared deletions within group media conditions
"""Finds the removed reactions for a given sample set specified by the accession number, number of samples, and control list 
(0 = control, 1 = experimental)"""
def findDeletedReactions(assent,CSmodel,numSamples,control,solver = 'gurobi',path= None):
    reactKOs = [set(), set()]
    for x in [z + 1 for z in range(numSamples)]:
        tempModel = CSmodel.copy()#create copy of model
        fileName = join(path, assent + '_' + str(x) + '.txt')
        statesFile = open(fileName, 'r')
        geneStates = statesFile.readlines()
        geneStates = [y.split(" ") for y in geneStates]#get gene states
        geneNames = []
        for y in geneStates:
            if int(y[1]) == 0: #find "off" genes
                geneNames.append(y[0])
        tempModel.solver = solver  # change solver
        if not (len(geneNames) == 0):
            cobra.manipulation.delete_model_genes(tempModel, geneNames)# create CS model
            temp = set()
            for r in tempModel.reactions: #find removed reactions
                if r.bounds == (0,0):
                    temp.add(r.id)
            if len(reactKOs[control[x-1]]) == 0:#add on removed reactions
               reactKOs[control[x-1]] = temp.copy()
            else:
               reactKOs[control[x-1]] |= temp.copy()#make union of sets

    return reactKOs

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
def createEssentialGeneDataAssent(assent,CSmodel,numSamples =1 ,path=None,solver='gurobi'):

    """For use in create CS model from geneStates.txt file in matlab"""
    newFile = open(('EssGene' + assent + '.txt'), 'w')
    for x in [x + 1 for x in range(numSamples)]:
        tempModel = CSmodel.copy()
        fileName = join(path, assent + '_' + str(x) + '.txt')
        statesFile = open(fileName,'r')
        geneStates = statesFile.readlines()
        geneStates = [y.split(" ") for y in geneStates]
        geneNames = []
        for y in geneStates:
            if int(y[1]) == 0:
                geneNames.append(y[0])
        tempModel.solver = solver  # load model
        if not(len(geneNames) == 0):
            cobra.manipulation.delete_model_genes(tempModel,geneNames)
        res = cobra.flux_analysis.single_gene_deletion(tempModel)  # perform single gene deletion simulation
        res.sort_index()  # sort the results of the simulation
        results = []
        for y in res.to_records(index=True):
            results.append(tuple(y))  # gather geneName,flux tuples
        newFile.write('\n-----------------------\n' + assent + '_' + str(x) + '\n')  # write sample label
        [newFile.write(str(x[0]) + ' ' + str(x[1]) + '\n') for x in results]  # write results to file
        statesFile.close()

    newFile.write('\n-----------------------\n')  # write final marker
    newFile.close()
"""Same as above but with reaction deletions"""
def createSingleReactionDeletionData(assent,model,numSamples=1,path=None,solver = 'gurobi'):
    newFile = open(('EssReact' + assent + '.txt'), 'w')
    for x in [x + 1 for x in range(numSamples)]:
        CSmodel = model.copy()
        fileName = join(path, assent + '_' + str(x) + '.txt')
        statesFile = open(fileName, 'r')
        geneStates = statesFile.readlines()
        geneStates = [y.split(" ") for y in geneStates]
        geneNames = []
        for y in geneStates:
            if int(y[1]) == 0:
                geneNames.append(y[0])
        CSmodel.solver = solver  # load model
        if not (len(geneNames) == 0):
            cobra.manipulation.delete_model_genes(CSmodel, geneNames)
        res = cobra.flux_analysis.single_reaction_deletion(CSmodel)  # perform single gene deletion simulation
        res.sort_index()  # sort the results of the simulation
        results = []
        for y in res.to_records(index=True):
            results.append(tuple(y))  # gather geneName,flux tuples
        newFile.write('\n-----------------------\n' + assent + '_' + str(x) + '\n')  # write sample label
        [newFile.write(str(x[0]) + ' ' + str(x[1]) + '\n') for x in results]  # write results to file
        statesFile.close()

    newFile.write('\n-----------------------\n')  # write final marker
    newFile.close()
"""Performs loopless FVA analysis of sample set. Uses the gene states given by accession#_sample#.txt to create CS model then perofrms analysis"""
def perform_fva(assent,model,numSamples=1,path = None,solver= 'gurobi'):
    totalFVARes = []
    for x in [x + 1 for x in range(numSamples)]:
        CSmodel = model.copy()
        fileName = join(path, assent + '_' + str(x) + '.txt')
        statesFile = open(fileName, 'r')
        geneStates = statesFile.readlines()
        geneStates = [y.split(" ") for y in geneStates]
        geneNames = []
        for y in geneStates:
            if int(y[1]) == 0:
                geneNames.append(y[0])
        CSmodel.solver = solver  # load model
        if not (len(geneNames) == 0):
            cobra.manipulation.delete_model_genes(CSmodel, geneNames)
        res = cobra.flux_analysis.flux_variability_analysis(CSmodel,loopless=True)  # perform single gene deletion simulation
          # gather geneName,flux tuples
        res.to_csv('FVA'+assent+'_'+str(x)+'.csv')
        statesFile.close()
    return totalFVARes

"""
Same as above for singel gene and reaction deletions but uses direct model input instead of creating within the function
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

def createEssentialReactionModel(model,label,solver= 'gurobi'):
    model.solver = solver
    res = cobra.flux_analysis.single_reaction_deletion(model)  # perform simulation
    newFile = open(('EssReact' + label + '.txt'), 'w')
    res.sort_index()
    results = []
    for x in res.to_records(index=True):
        results.append(tuple(x))  # gather geneName flux level tuple
    [newFile.write('\n' + str(x[0]) + ' ' + str(x[1])) for x in results]  # write results
    newFile.write('\n')
    newFile.close()



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
        if len(self.differences) > 1:
           # print (len(self.differences))#TODO add in null data to the differences array to show more accurate statistic; multiply by control coeff pad with 0
            _,pval = ttest_1samp(self.differences,0.0)
        else:
            pval = 1.0
        return self.name, self.count ,mean(self.differences),std(self.differences),sum([abs(x) for x in self.differences]),pval
    def areEqual(self,x):
        if x.name == self.name:
            return 1
        else:
            return 0
    """Combines two GeneHits objects"""
    def combine(self,x):
        self.count += x.count
        self.differences = self.differences + x.differences

"""
Class for the characterization of essential genes for a general metabolic model. Processes the data in the EssGene file. 
Data is acessible through the getData() function where the genes and corresponding fluxes are stored as a dictionary. For
use in comparision to CS models
"""
class ComparisionGene:
    """creates object that contains the results of a single gene or reaction deletion. Uses the label to find .txt file"""
    def __init__(self,label,information,type = 'Gene'):
        self.normalFile = open('Ess'+type+label+'.txt', 'r')
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
. Allows for comparision of data from the file with a general model data. Can handle the control and experimental data as well as media conditons

"""
class CSGene:
    def __init__(self,accession,control,media =1,information='NO INFORMATION',type = 'Gene'):
        self.accession = accession
        self.information = information
        self.file = open('Ess'+type+self.accession+'.txt','r')
        self.dataTot = self.file.readlines() #collect values in list of strings
        self.control = control #define array of control information e.g. [1,0,0,0,1,1]
        self.media = media #define media value
        try:
            while(True):
                self.dataTot.remove('\n') #remove lines
        except:
            pass
        for x in range(len(self.dataTot)):
            self.dataTot[x] = self.dataTot[x][:-1] #remove eol characters
        self.file.close()
    """Splits the name and flux, and separates based on experimental or control conditions"""
    def processSamples(self):
        self.sampleEG = [] #2d list for holding gene, flux
        tempData = []
        for x in self.dataTot:
            if x[0] == '-': #if starting new sample
                if len(self.sampleEG) == 0:
                    self.sampleEG = [tempData] #if on the second sample add all to the sampleEG list
                else:
                    tempData2 = [tempData[0]]
                    for y in [x+1 for x in range(len(tempData)-1)]:
                        spliter = tempData[y].split(" ") #split values
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
    """find differences in name flux pairs to some other group (generalGene) and does so only with the group and media specified by group (1,0) or media (0,1,2).
    If Media = 3 then it works irrespective of media conditions"""
    def findChangedFluxGenes(self,generalGene,group,media = 3):
        results = []
        data = self.sampleEG
        i = 0
        for x in data:
            coeff = self.control[i]
            for y in data[x]:
                if not(data[x][y] == generalGene[y]) and coeff == group and (self.media == media or media == 3) :#make sure sample matches
                    if len(results) == 0:
                        results.append(GeneHits(data[x][y] - generalGene[y],y))
                    else:
                        count = 0
                        for z in results:
                            if z.name == y:
                                z.combine(GeneHits((data[x][y] - generalGene[y]),y))#merge geneHist object
                                count +=1
                        if count == 0:
                            results.append(GeneHits((data[x][y] - generalGene[y]),y))
            i +=1
        numE = 0#pad results with 0 to represent 0 flux change in samples
        if self.media == media or media == 3:
            for x in self.control:
                if x == group:
                    numE += 1
            for x in range(len(results)):
                if len(results[x].differences) > numE:
                    temp = array(results[x].differences)
                    temp.sort()
                    temp = list(temp)
                    results[x].differences = temp[0:numE-1]
                while len(results[x].differences) < numE:
                    results[x].differences.append(0.0)

        return results,numE











