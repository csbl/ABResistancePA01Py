from EssentialGeneFinder import *
from copy import deepcopy


def creationOfEssentialGeneData():
    assent = ['GSE90620']#,'GDS4244','GDS3572','GSE30021','GSE65870']
    nsamples = [12]#,12,6,9,6]

    for i in range(len(assent)):

        createEssentialGeneDataAssent(assent[i],'.mat',nsamples[i],'C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')
      #  model = cobra.io.read_sbml_model("iPAE1146.xml")
       # createEssentialGeneModel(model,'IPAE1146')


creationOfEssentialGeneData()

model = cobra.io.read_sbml_model("iPAE1146.xml")
createEssentialGeneModel(model,'IPAE1146')

normal = ComparisionGene('IPAE1146','This is a test')


reference = {'GSE90620':[1,1,1,0,0,0,0,0,0,0,0,0],'GDS3572':[0,0,0,1,1,1],'GSE65870':[0,0,0,0,0,0],'GDS4244':[1,1,1,1,1,1,0,0,0,0,0,0],'GSE30021':[1,1,1,1,1,1,0,0,0]}

CSD = []
"""

for x in reference:
    temp = CSGene(x,reference[x])
    temp.processSamples()
    CSD.append(temp)

uniqueData = list()

for x in CSD:
    temp,c,b = x.findUniqueFrom(normal.getData())
    uniqueData.append(temp)
    pings = list()
    count = 0

pings = findNumberofHits(uniqueData)

[j.letsPrint() for j in pings]
"""
