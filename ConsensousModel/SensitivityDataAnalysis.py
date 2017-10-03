import pandas

#Load in gene lists

geneListF = open('RUniqGenes.txt', 'r')
geneListR = [x[:-1] for x in geneListF.readlines()]
geneListF.close()
geneListF = open('CUniqGenes.txt', 'r')
geneListC = [x[:-1] for x in geneListF.readlines()]
geneListF.close()


results = dict()
orderOfDel = [[[x],[]] for x in geneListC] + [[[],[x]] for x in geneListR]#delete impactful genes sequentially

data = pandas.DataFrame.from_csv("sensitivityResultsFullData.csv")
data.columns = ['KOReactions','FVA','GeneDeletion','RxnDeletion','SensitiveGene','ResistantGene']

# Key for data frame column lists: 0:baseline[s] - experimental[s], 1: experimental[s] - baseline[s], 2:baseline[r] - experimental[r], 1: experimental[r] - baseline[r]



print data