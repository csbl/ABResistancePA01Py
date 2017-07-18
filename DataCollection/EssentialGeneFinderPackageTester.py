from EssentialGeneFinder import *
import scipy.stats
from pandas import DataFrame,Series


def creationOfEssentialGeneData():
    assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
    nsamples = [12,12,6,9,6]
    model = cobra.io.read_sbml_model("iPAE1146.xml")
    createEssentialGeneModel(model,'IPAE1146')

    for i in range(len(assent)):
        createEssentialGeneDataAssent(assent[i],model,nsamples[i],'C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')





#creationOfEssentialGeneData()

normal = ComparisionGene('IPAE1146','This is a test')

reference = {'GSE90620':[0,0,0,1,1,1,1,1,1,1,1,1],'GDS3572':[1,1,1,0,0,0],'GSE65870':[1,1,1,1,1,1],'GDS4244':[0,0,0,0,0,0,1,1,1,1,1,1],'GSE30021':[0,0,0,0,0,0,1,1,1]}
referenceMed =  {'GSE90620':1,'GDS3572':0,'GSE65870':2,'GDS4244':1,'GSE30021':1}
# PBM = 2
# LB = 1
# MH = 0
CSD = []

for x in reference:
    temp = CSGene(x,reference[x],referenceMed[x])
    temp.processSamples()
    CSD.append(temp)

results = []
i=0
sampleCount = [0,0]
typeOfSample = [0,1]
controlSample = DataFrame(columns = ['Name', 'Count', 'Mean', 'Std', 'Sum', 'Pvalue', 'GEM_FBM', 'CS_FBM'])
experimentalSample = DataFrame(columns = ['Name', 'Count', 'Mean', 'Std', 'Sum', 'Pvalue', 'GEM_FBM', 'CS_FBM'])
controlMedia = DataFrame()
experimentalMedia = DataFrame()
data = [controlSample,experimentalSample]
mediaData = [controlMedia,experimentalMedia]
media = ['MH','LB','PBM']
for g in typeOfSample:
    for m in [0,1,2,3]:
        results = []
        i = 0
        j = 0
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
            sampleCount[g] += tempCount

        types = {1: 'Experimental', 0: 'Control'}
        if not m == 3:
            tempMediaDataFrame = DataFrame(columns = ['Name',media[m]+' Count',media[m] + ' Sum'])
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
        y.append(t_testUnpaired_fromSum(m,sd**2,sampleCount[1],con["Mean"].values,con['Std'].values**2,sampleCount[0]))
    except:
        y.append(0.0)
data[1].insert(len(data[1].columns.values),"PVal(cont)", y)
i = 0
for x in data:
    print x
    print sampleCount[i]
    x.to_csv(types[i]+'.txt',sep = " ",float_format="%.5f",index=False)
    i+=1


