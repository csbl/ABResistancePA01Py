from ABResistanceAnalysis import *
from pandas import DataFrame,Series
from matplotlib import pyplot as plt
import numpy
import matplotlib


"""Script for creation of flux variability analysis data and subsequent comparision between control and experimental groups
"""
def creationOfFVAData():
    assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
    nsamples = [12,12,6,9,6]
    model = cobra.io.read_sbml_model("iPAE1146.xml")
    createEssentialReactionModel(model,'IPAE1146')
    comp_res = cobra.flux_analysis.flux_variability_analysis(model)
    comp_res.to_csv('FVAIPAE1146.csv')
    for i in range(len(assent)):
        perform_fva(assent[i],model,nsamples[i],'C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')




#creationOfFVAData()#uncomment if FVA data needs to be created for current models
model = cobra.io.read_sbml_model("iPAE1146.xml")#load general model
assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
reference = {'GSE90620':[0,0,0,1,1,1,1,1,1,1,1,1],'GDS3572':[1,1,1,0,0,0],'GSE65870':[1,1,1,1,1,1],'GDS4244':[0,0,0,0,0,0,1,1,1,1,1,1],'GSE30021':[0,0,0,0,0,0,1,1,1]}
nsamples = [12,12,6,9,6]
mini = 0#set bounds for matrix correspond to a bounding limit of -1000 and 1000
maxi = 2000
n = 5 #number of different sample sets
compModelFVA =DataFrame.from_csv('FVA'+'IPAE1146'+'.csv')
types = {1: 'Experimental FVA', 0: 'Control FVA'}
xbins = len(compModelFVA.index.values) #find number of reactions
lb = abs(min([x.lower_bound for x in model.reactions]))#find lower bound
ub = abs(max([x.upper_bound for x in model.reactions]))#find upper bound
ysize = lb + ub #create sizing
ybins = maxi-mini
yblim = 900 #shown limit on graph (-100)
ytlim = 2000#shown upper limit (1000)
geneMap = Series(range(xbins),compModelFVA.index.values).to_dict()#map reactions to numbers
geneMapInv = Series(compModelFVA.index.values,range(xbins)).to_dict()#inverse mapping
topHits = dict()
FVARes = []
for group in [0,1]:#iterate through groups
    #fig = plt.figure(figsize = (4,6),dpi = 100)
    resultMat = numpy.zeros((xbins,ysize))#matrix for results
    for i in range(len(assent)):
        for j in [x+1 for x in range(nsamples[i])]:
            if reference[assent[i]][j-1] == group:
                data = DataFrame.from_csv('FVA'+assent[i]+'_'+str(j)+'.csv')
                for index,row in data.iterrows():#iterate through FVA data for sample
                    lb2 = row['minimum']
                    ub2 = row['maximum']
                    if lb2 == 0 and ub2 == 0:
                        resultMat[geneMap[index],lb] += 1 #if rxn is removed
                    span = range(int(ub2)) + lb2 #find range of flux values
                    for y in span:#bin the range
                        #if y == 0 and len(span) == 1:
                        resultMat[geneMap[index],y+lb] += 1
    topHits[types[group]] = numpy.greater(resultMat,.7*numpy.max(numpy.max(resultMat))).astype(int)#find the range common to >70% of samples in a group
    resHits = topHits[types[group]][mini:maxi]#pull out data
    plt.matshow(resHits.transpose(),interpolation = 'none',cmap = 'binary')#plot result
    #plt.colorbar(fraction=0.046, pad=0.04)
    plt.xticks(rotation = 45)
    #plt.xticks(range(xbins),[geneMapInv[x] for x in range(xbins)])
    plt.xticks()
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        top='off',  # ticks along the top edge are off
        labeltop='off')  # labels along the bottom edge are off
    plt.yticks(numpy.linspace(mini,maxi,n),[x-lb for x in numpy.linspace(mini,maxi,n)])
    plt.ylim((yblim,1600))
    plt.xlabel('Reaction')
    plt.ylabel('Flux')
    plt.title(types[group]+' '+'Top Hits', y = 1.15)
    plt.savefig(types[group]+'topHits'+'.png',bbox_inches='tight')
    resPlot = resultMat[:][mini:maxi]#plot total differences
    plt.matshow(resPlot.transpose(),interpolation = 'none',cmap = 'binary')
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
diff = abs(FVARes[0] - FVARes[1])#find difference in matrices for two groups to quantify FVA differences
print numpy.shape(diff)
diffSum =  sum(diff,axis = 0)
change = Series(numpy.abs(diffSum),[geneMapInv[x] for x in range(xbins)])
change.to_csv('FVADiff.csv') #output change

change = Series(sum(numpy.abs(numpy.subtract(topHits[types[0]],topHits[types[1]])),axis= 1),[geneMapInv[x] for x in range(xbins)])
change.to_csv('FVADiffTopHits.csv') #find difference in > 70% range
plt.show()
change = DataFrame(sum(numpy.abs(numpy.subtract(topHits[types[0]],topHits[types[1]])),axis= 1),index = [geneMapInv[x] for x in range(xbins)], columns = ['DiffFlux'])
change = change.query('DiffFlux > 0') #tabulate top differences
print change
rxns = [geneMap[x] for x in change.index.values]#convert to rxn names
lbounds = [[],[]]
ubounds = [[],[]]

## find corresponding bounds
for y in [0, 1]:
    for x in rxns:
        flag = 0
        for z in range(2000):
            if( not (topHits[types[y]][x,z] == 0) and flag == 0):#grab lowerbound
                lbounds[y].append(z- lb)
                flag += 1
            if (topHits[types[y]][x,z] == 0 and flag == 1):#grab upper bound
                ubounds[y].append(z-lb - 1)
                flag += 1
            if z == 1999 and flag == 1:#if going to end grab upper bound
                ubounds[y].append(1000)
            if z == 1999 and flag == 0:#if never starts grab lower bound
                lbounds[y].append(0)
                ubounds[y].append(0)
#output dataframe to csv
change.insert(len(change.columns.values),'lboundC',value = Series(lbounds[0],index = [geneMapInv[x] for x in rxns]))
change.insert(len(change.columns.values),'uboundC',Series(ubounds[0],index = [geneMapInv[x] for x in rxns]))
change.insert(len(change.columns.values),'lboundE',Series(lbounds[1],index = [geneMapInv[x] for x in rxns]))
change.insert(len(change.columns.values),'uboundE',Series(ubounds[1],index = [geneMapInv[x] for x in rxns]))
change.to_csv('FVATopHitsWFluxes.csv')




