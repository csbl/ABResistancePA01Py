from EssentialGeneFinder import *
from pandas import DataFrame
import cobra
import cobra.manipulation

def creationOfDeletedReaction():
    assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
    reference = {'GSE90620': [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1], 'GDS3572': [1, 1, 1, 0, 0, 0],
                 'GSE65870': [1, 1, 1, 1, 1, 1], 'GDS4244': [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                 'GSE30021': [0, 0, 0, 0, 0, 0, 1, 1, 1]}
    referenceMed = {'GSE90620': 1, 'GDS3572': 0, 'GSE65870': 2, 'GDS4244': 1, 'GSE30021': 1}
    nsamples = [12,12,6,9,6]
    deletions = [[set(),set()],[set(),set(),set()]]
    model = cobra.io.read_sbml_model("iPAE1146.xml")
    createEssentialReactionModel(model,'IPAE1146')
    for i in range(len(assent)):
        res = findDeletedReactions(assent[i],model,nsamples[i],reference[assent[i]],'gurobi','C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')
        if not len(res[0]) == 0 : deletions[0][referenceMed[assent[i]]] = deletions[0][referenceMed[assent[i]]].union(res[0])
        if not len(res[1]) == 0 : deletions[1][referenceMed[assent[i]]] = deletions[1][referenceMed[assent[i]]].union(res[1])
        if referenceMed[assent[i]] == 1:
            print res[1]
            print  " \n"
       # print deletions[1][1]
    #for x in deletions[1]:
      #  print x
       # print "\n\n"
    return [deletions[0][1] & deletions[0][1],deletions[1][0] | deletions[1][1] | deletions[1][2]]
def findDeletedReactions(assent,CSmodel,numSamples,control,solver = 'gurobi',path= None):
    newFile = open(('EssGene' + assent + '.txt'), 'w')
    for x in [x + 1 for x in range(numSamples)]:
        tempModel = CSmodel.copy()
        fileName = join(path, assent + '_' + str(x) + '.txt')
        statesFile = open(fileName, 'r')
        geneStates = statesFile.readlines()
        geneStates = [y.split(" ") for y in geneStates]
        geneNames = []
        for y in geneStates:
            if int(y[1]) == 0:
                geneNames.append(y[0])
        tempModel.solver = solver  # load model
        reactKOs = [set(),set()]
        if not (len(geneNames) == 0):
            cobra.manipulation.delete_model_genes(tempModel, geneNames)
            temp = set()
            for r in tempModel.reactions:
                if r.bounds == (0,0):
                    temp.add(r.id)
            if len(reactKOs[control[x-1]]) == 0:
               reactKOs[control[x-1]] = temp
            else:
                reactKOs[control[x-1]] = reactKOs[control[x-1]] | temp
        return reactKOs


print creationOfDeletedReaction()
"""
normal = ComparisionGene('IPAE1146','This is a test',type = 'React')

reference = {'GSE90620':[0,0,0,1,1,1,1,1,1,1,1,1],'GDS3572':[1,1,1,0,0,0],'GSE65870':[1,1,1,1,1,1],'GDS4244':[0,0,0,0,0,0,1,1,1,1,1,1],'GSE30021':[0,0,0,0,0,0,1,1,1]}
referenceMed =  {'GSE90620':1,'GDS3572':0,'GSE65870':2,'GDS4244':1,'GSE30021':1}
# PBM = 2
# LB = 1
# MH = 0
CSD = []
for x in reference:
    temp = CSGene(x,reference[x],referenceMed[x],type = 'React')
    temp.processSamples()
    CSD.append(temp)


controlSample = DataFrame(columns = ['Name', 'Count', 'Mean', 'Std', 'Sum', 'Pvalue', 'GEM_FBM', 'CS_FBM'])
experimentalSample = DataFrame(columns = ['Name', 'Count', 'Mean', 'Std', 'Sum', 'Pvalue', 'GEM_FBM', 'CS_FBM'])
data = [controlSample,experimentalSample]
controlMedia = DataFrame()
experimentalMedia = DataFrame()
mediaData = [controlMedia,experimentalMedia]
media = ['MH','LB','PBM']
sampleCount = [0,0]
typeOfSample = [0,1]
results = []
i=0
for g in typeOfSample:
    for m in [0,1,2,3]:
        results = []
        i = 0
        j=0
        sampleCount[g] = 0
        for x in CSD:
            temp,tempCount = x.findChangedFluxGenes(normal.getData(),g,m)
            s = len(results)
            if i==0:
                results = temp[:]
            else:
                for y in temp:
                    count = 0
                    for z in range(s):
                        if results[z].areEqual(y):
                            results[z].combine(y)
                            count += 1
                            break
                    if count == 0:
                        results = results + [y]
            i+=1
            sampleCount[g]+= tempCount

        types = {1: 'ExperimentalReact', 0: 'ControlReact'}
        if not m == 3:
            tempMediaDataFrame = DataFrame(columns=['Name', media[m] + ' Count', media[m] + ' Sum'])
        for x in results:
            name,count,mean,std,sum,pval = x.letsPrint()
            if m == 3:
                if( abs(mean) > .001 ):
                    data[g].loc[j] = [name,count,mean,std,sum,pval,str(normal.getData()[x.name]),normal.getData()[x.name]+ mean]
                    j += 1
            else:
                if( abs(mean) > .001 ):
                    tempMediaDataFrame.loc[j] = [name,count,sum]
                    j += 1
        if m == 3:
            data[g] = data[g].merge(mediaData[g],on = 'Name')
        elif m==0:
            mediaData[g] = tempMediaDataFrame.copy()
        elif not(tempMediaDataFrame.empty):
            mediaData[g] = mediaData[g].merge(tempMediaDataFrame,on = 'Name')

y = []
for z,m,sd in zip(data[1]["Name"].values,data[1]["Mean"].values,data[1]["Std"].values):
    try:
        con = controlSample.loc[controlSample['Name'] == z]
        y.append(t_testUnpaired_fromSum(m,sd**2,sampleCount[1],con['Mean'].values,con['Std'].values**2,sampleCount[0]))
    except:
        y.append(0.0)
data[1].insert(len(data[1].columns.values),"PVal(cont)", y)
i = 0
for x in data:
    print x
    print sampleCount[i]
    x.to_csv(types[i]+'.txt',sep = " ",float_format="%.5f",index=False)
    i+=1
"""



