from EssentialGeneFinder import *
from pandas import DataFrame,Series
from matplotlib import pyplot as plt
import numpy
import matplotlib



def creationOfFVAData():
    assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
    nsamples = [12,12,6,9,6]
    model = cobra.io.read_sbml_model("iPAE1146.xml")
    createEssentialReactionModel(model,'IPAE1146')
    comp_res = cobra.flux_analysis.flux_variability_analysis(model)
    comp_res.to_csv('FVAIPAE1146.csv')
    for i in range(len(assent)):
        perform_fva(assent[i],model,nsamples[i],'C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')




#creationOfFVAData()

assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
reference = {'GSE90620':[0,0,0,1,1,1,1,1,1,1,1,1],'GDS3572':[1,1,1,0,0,0],'GSE65870':[1,1,1,1,1,1],'GDS4244':[0,0,0,0,0,0,1,1,1,1,1,1],'GSE30021':[0,0,0,0,0,0,1,1,1]}
nsamples = [12,12,6,9,6]
compModelFVA =DataFrame.from_csv('FVA'+'IPAE1146'+'.csv')
types = {1: 'ExperimentalFVA', 0: 'ControlFVA'}
xbins = len(compModelFVA.index.values)
ybins = 2000
geneMap = Series(range(xbins),compModelFVA.index.values).to_dict()
geneMapInv = Series(compModelFVA.index.values,range(xbins)).to_dict()
min = 0
max = 2000
n = 5
FVARes = []
for group in [0,1]:
    #fig = plt.figure(figsize = (4,6),dpi = 100)
    resultMat = numpy.zeros((xbins,ybins))
    for i in range(len(assent)):
        for j in [x+1 for x in range(nsamples[i])]:
            if reference[assent[i]][j-1] == group:
                data = DataFrame.from_csv('FVA'+assent[i]+'_'+str(j)+'.csv')
                for index,row in data.iterrows():
                    lb = row['minimum']
                    ub = row['maximum']
                    span = range(int(ub)) + lb
                    for y in span:
                        resultMat[geneMap[index],y+1000] += 1
    plt.matshow(resultMat.transpose()[min:max][:],interpolation = 'none',cmap = 'hot')
    plt.colorbar()
    plt.xticks(rotation = 45)
    #plt.xticks(range(xbins),[geneMapInv[x] for x in range(xbins)])
    plt.yticks(numpy.linspace(min,max,n),[x-1000 for x in numpy.linspace(min,max,n)])
    plt.xlabel('Gene')
    plt.ylabel('Flux')
    plt.title(types[group], y = 1.25)
    FVARes.append(resultMat.transpose().copy())

diff = FVARes[0] - FVARes[1]
diffSum =  sum(diff,axis = 0)
change = Series(numpy.abs(diffSum),[geneMapInv[x] for x in range(xbins)])
change.to_csv('FVADiff.csv')

plt.show()

"""
normal = ComparisionGene('IPAE1146','This is a test',type = 'React')

reference = {'GSE90620':[0,0,0,1,1,1,1,1,1,1,1,1],'GDS3572':[1,1,1,0,0,0],'GSE65870':[1,1,1,1,1,1],'GDS4244':[0,0,0,0,0,0,1,1,1,1,1,1],'GSE30021':[0,0,0,0,0,0,1,1,1]}

CSD = []

controlSample = DataFrame(columns = ['Name', 'Count', 'Mean', 'Std', 'Sum', 'Pvalue', 'GEM_FBM', 'CS_FBM'])
experimentalSample = DataFrame(columns = ['Name', 'Count', 'Mean', 'Std', 'Sum', 'Pvalue', 'GEM_FBM', 'CS_FBM'])
data = [controlSample,experimentalSample]
for x in reference:
    temp = CSGene(x,reference[x],type = 'React')
    temp.processSamples()
    CSD.append(temp)

results = []
i=0
sampleCount = [0,0]
typeOfSample = [0,1]
for g in typeOfSample:
    results = []
    i = 0
    j=0
    sampleCount[g] = 0
    for x in CSD:
        temp,tempCount = x.findChangedFluxGenes(normal.getData(),g)
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
    for x in results:
        name,count,mean,std,sum,pval = x.letsPrint()
        if( abs(mean) > .001 ):
            data[g].loc[j] = [name,count,mean,std,sum,pval,str(normal.getData()[x.name]),normal.getData()[x.name]+ mean]
            j+=1
    print data[g]
    footer = 'Total Number of %s Samples = %d' % (types[g],sampleCount[g])
    print footer
y = []
for z,m,sd in zip(data[1]["Name"].values,data[1]["Mean"].values,data[1]["Std"].values):
    try:
        con = controlSample.loc[controlSample['Name'] == z]
        y.append(t_testUnpaired_fromSum(m,sd**2,sampleCount[1],con['Mean'].values,con['Std'].values**2,sampleCount[0]))
    except:
        y.append(0.0)
    print y
data[1].insert(len(data[1].columns.values),"PVal(cont)", y)
i = 0;
for x in data:
    x.to_csv(types[i]+'.txt',sep = " ",float_format="%.5f",index=False)
    i+=1
"""