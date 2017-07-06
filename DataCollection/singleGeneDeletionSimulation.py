from EssentialGeneFinder import *
import scipy.stats
from pandas import DataFrame


def creationOfEssentialGeneData():
    assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
    nsamples = [12,12,6,9,6]
    model = cobra.io.read_sbml_model("iPAE1146.xml")
    createEssentialGeneModel(model,'IPAE1146')

    for i in range(len(assent)):
        createEssentialGeneDataAssent(assent[i],model,nsamples[i],'C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')
def t_testUnpaired_fromSum(y1,s1,n1,y2,s2=-1,n2= -1,a = .05):
    if s2 == -1: s2 = s1
    if n2 == -1: n2 = n1
    T = (y1-y2)/sqrt(s1/n1 + s2 / n2)
    v = ((s1/n1+s2/n2)**2)/(((s1/n1)**2)/(n1-1)+((s2/n2)**2)/(n2-1))
    return scipy.stats.t.sf(T,v)*2




#creationOfEssentialGeneData()

normal = ComparisionGene('IPAE1146','This is a test')

reference = {'GSE90620':[0,0,0,1,1,1,1,1,1,1,1,1],'GDS3572':[1,1,1,0,0,0],'GSE65870':[1,1,1,1,1,1],'GDS4244':[0,0,0,0,0,0,1,1,1,1,1,1],'GSE30021':[0,0,0,0,0,0,1,1,1]}

CSD = []

for x in reference:
    temp = CSGene(x,reference[x])
    temp.processSamples()
    CSD.append(temp)

results = []
i=0
sampleCount = 0
typeOfSample = [0,1]
controlSample = DataFrame(columns = ['Name', 'Count', 'Mean', 'Std', 'Sum', 'Pvalue', 'GEM_FBM', 'CS_FBM'])
experimentalSample = DataFrame(columns = ['Name', 'Count', 'Mean', 'Std', 'Sum', 'Pvalue', 'GEM_FBM', 'CS_FBM'])
data = [controlSample,experimentalSample]
for g in typeOfSample:
    results = []
    i = 0
    j = 0
    sampleCount = 0
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
        sampleCount += tempCount

    types = {1: 'Experimental', 0: 'Control'}
    file = open(types[g]+'.txt','w')
    headers = 'Name Count Mean Std Sum Pvalue GEM_FBM CS_FBM'
    print headers
    file.write(headers+'\n')
    for x in results:
        name,count,mean,std,sum,pval = x.letsPrint()
        if( abs(mean) > .001 ):
            resulting = '%s  %d  %.12f  %.12f  %.12f  %.12f  %s  %.4f' % (name,count,mean,std,sum,pval,str(normal.getData()[x.name]),normal.getData()[x.name]+ mean)
            file.write(resulting+'\n')
            data[g].loc[j] = [name,count,mean,std,sum,pval,str(normal.getData()[x.name]),normal.getData()[x.name]+ mean]
            print resulting
            j += 1

    footer = 'Total Number of %s Samples = %d' % (types[g],sampleCount)
    print footer
    file.close()

file = open('ExperimentalWithTTest.txt','w')
file.write(headers+' PVal (ctc)'+'\n')
y = experimentalSample.values.tolist()
for z in y:
    resulting = '%s  %d  %.12f  %.12f  %.12f  %.12f  %s  %.4f' % tuple(z)
    try:
        con = controlSample.loc[controlSample['Name'] == z[0]]
        file.write(resulting + '    %.4f' % t_testUnpaired_fromSum(z[2],z[3]**2,sampleCount,con['Mean'].values,con['Std'].values,18) + '\n')
    except:
        file.write(resulting + '    0.0' + '\n')
file.close()

