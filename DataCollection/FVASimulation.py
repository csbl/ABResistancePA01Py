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




creationOfFVAData()
model = cobra.io.read_sbml_model("iPAE1146.xml")
assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
reference = {'GSE90620':[0,0,0,1,1,1,1,1,1,1,1,1],'GDS3572':[1,1,1,0,0,0],'GSE65870':[1,1,1,1,1,1],'GDS4244':[0,0,0,0,0,0,1,1,1,1,1,1],'GSE30021':[0,0,0,0,0,0,1,1,1]}
nsamples = [12,12,6,9,6]
mini = 0
maxi = 2000
n = 5
compModelFVA =DataFrame.from_csv('FVA'+'IPAE1146'+'.csv')
types = {1: 'ExperimentalFVA', 0: 'ControlFVA'}
xbins = len(compModelFVA.index.values)
lb = abs(min([x.lower_bound for x in model.reactions]))
ub = abs(max([x.upper_bound for x in model.reactions]))
ysize = lb + ub
ybins = maxi-mini
yblim = 900
ytlim = 2000
geneMap = Series(range(xbins),compModelFVA.index.values).to_dict()
geneMapInv = Series(compModelFVA.index.values,range(xbins)).to_dict()
topHits = dict()
FVARes = []
for group in [0,1]:
    #fig = plt.figure(figsize = (4,6),dpi = 100)
    resultMat = numpy.zeros((xbins,ysize))
    for i in range(len(assent)):
        for j in [x+1 for x in range(nsamples[i])]:
            if reference[assent[i]][j-1] == group:
                data = DataFrame.from_csv('FVA'+assent[i]+'_'+str(j)+'.csv')
                for index,row in data.iterrows():
                    lb2 = row['minimum']
                    ub2 = row['maximum']
                    span = range(int(ub2)) + lb2
                    for y in span:
                        resultMat[geneMap[index],y+lb] += 1
    topHits[types[group]] = numpy.equal(resultMat,numpy.max(numpy.max(resultMat))).astype(int)
    resHits = topHits[types[group]][mini:maxi]
    plt.matshow(resHits.transpose(),interpolation = 'none',cmap = 'hot')
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.xticks(rotation = 45)
    #plt.xticks(range(xbins),[geneMapInv[x] for x in range(xbins)])
    plt.yticks(numpy.linspace(mini,maxi,n),[x-lb for x in numpy.linspace(mini,maxi,n)])
    plt.ylim((yblim,ytlim))
    plt.xlabel('Reaction')
    plt.ylabel('Flux')
    plt.title(types[group]+'Top Hits', y = 1.15)
    resPlot = resultMat[:][mini:maxi]
    plt.matshow(resPlot.transpose(),interpolation = 'none',cmap = 'hot')
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.xticks(rotation = 45)
    #plt.xticks(range(xbins),[geneMapInv[x] for x in range(xbins)])
    plt.yticks(numpy.linspace(mini,maxi,n),[x-lb for x in numpy.linspace(mini,maxi,n)])
    plt.ylim((yblim,ytlim))
    plt.xlabel('Reaction')
    plt.ylabel('Flux')
    plt.title(types[group], y = 1.15)
    FVARes.append(resultMat.transpose().copy())
    plt.savefig(types[group]+'.png',bbox_inches='tight')
diff = abs(FVARes[0] - FVARes[1])
print numpy.shape(diff)
diffSum =  sum(diff,axis = 0)
change = Series(numpy.abs(diffSum),[geneMapInv[x] for x in range(xbins)])
change.to_csv('FVADiff.csv')

change = Series(sum(numpy.abs(numpy.subtract(topHits[types[0]],topHits[types[1]])),axis= 1),[geneMapInv[x] for x in range(xbins)])
change.to_csv('FVADiffTopHits.csv')
plt.show()
change = DataFrame(sum(numpy.abs(numpy.subtract(topHits[types[0]],topHits[types[1]])),axis= 1),index = [geneMapInv[x] for x in range(xbins)], columns = ['DiffFlux'])
change = change.query('DiffFlux > 0')
print change
rxns = [geneMap[x] for x in change.index.values]
lbounds = [[],[]]
ubounds = [[],[]]


for y in [0, 1]:
    for x in rxns:
        flag = 0
        for z in range(2000):
            if( not (topHits[types[y]][x,z] == 0) and flag == 0):
                lbounds[y].append(z- lb)
                flag += 1
            if (topHits[types[y]][x,z] == 0 and flag == 1):
                ubounds[y].append(z-lb - 1)
                flag += 1
            if z == 1999 and flag == 1:
                ubounds[y].append(1000)
            if z == 1999 and flag == 0:
                lbounds[y].append(0)
                ubounds[y].append(0)

change.insert(len(change.columns.values),'lboundC',value = Series(lbounds[0],index = [geneMapInv[x] for x in rxns]))
change.insert(len(change.columns.values),'uboundC',Series(ubounds[0],index = [geneMapInv[x] for x in rxns]))
change.insert(len(change.columns.values),'lboundE',Series(lbounds[1],index = [geneMapInv[x] for x in rxns]))
change.insert(len(change.columns.values),'uboundE',Series(ubounds[1],index = [geneMapInv[x] for x in rxns]))
change.to_csv('FVATopHitsWFluxes.csv')




