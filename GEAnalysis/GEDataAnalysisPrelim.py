import numpy
from DataCollection import EssentialGeneFinder
import pandas
class geneHits:
    def __init__(self , name):
        self.name = name
        self.count = 1
       # self.index = name[:6]
    def addToCount(self):
        self.count = self.count + 1
    def letsPring(self):
        print self.name + ': %d' % self.count



accension = ['GDS3572','GDS4244','GSE30021','GSE65870','GSE90620']
normalFile = open('EssGeneIPAE1146.txt','r')

normalData = normalFile.readlines()

for x in normalData:
    if x == '\n' or x.startswith('-') == 1:
        normalData.remove(x)

data = []
tempData = []
for x in accension:
    file = open('EssGene'+x+'_.txt','r')
    #print 'EssGene' + x + '_.txt'
    dataTot = file.readlines()
    #print dataTot
    for x in dataTot:
        if x == '\n' or x.startswith('G')==1:
            dataTot.remove(x)
        elif x.startswith('-') == 1:
            data.append(tempData)
            tempData = []
        else:
            tempData.append(x[:-1])


uniqueData = []
for x in data:
    tempData = []
    for g in x:
        if EssentialGeneFinder.binary_search(normalData,g) == -1:
            tempData.append(g)
    uniqueData.append(tempData)

pings = []
count = 0

for x in uniqueData:
    for g in x:
        for r in x[g]:
            count = 0
            if len(pings) == 0:
                pings = [geneHits(g)]
            else:
                for j in pings:
                    if j.name == g:
                        count = count + 1
                if count ==0:
                    pings.append(geneHits(g))
                else:
                    for j in pings:
                        if j.name == g:
                            j.addToCount()

[j.letsPring() for j in pings]




