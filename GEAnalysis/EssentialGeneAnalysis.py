from DataCollection import EssentialGeneFinder

class GeneHits:
    def __init__(self , name):
        self.name = name
        self.count = 1
       # self.index = name[:6]
    def addToCount(self):
        self.count = self.count + 1
    def letsPring(self):
        print self.name + ': %d' % self.count
class ComparisionModel:
    def __init__(self,label,information):
        self.normalFile = open('EssGene'+label+'.txt', 'r')
        self.normalData = self.normalFile.readlines()
        for x in self.normalData:
            if x == '\n' or x.startswith('-') == 1:
                self.normalData.remove(x)
        self.normalFile.close()
        self.information = information
    def getData(self):
        return self.normalData
class