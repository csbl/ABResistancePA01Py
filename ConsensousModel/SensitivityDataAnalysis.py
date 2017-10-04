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

data = pandas.DataFrame.from_csv("sensitivityResultsFulData.csv")
data.columns = ['KOReactions','FVA','GeneDeletion','RxnDeletion','SensitiveGene','ResistantGene']

# Key for data frame column lists: 0:baseline[s] - experimental[s], 1: experimental[s] - baseline[s], 2:baseline[r] - experimental[r], 1: experimental[r] - baseline[r]
total = pandas.DataFrame(index = range(15))

for column in data.columns.values[:-2]:
    test = data[column].tolist()
   # print eval(test[0])
    bes = list()
    ebs = list()
    ber = list()
    ebr = list()

    for x in test:
        x = eval(x)
        bes.append(len(x[0]))
        ebs.append(len(x[1]))
        ber.append(len(x[2]))
        ebr.append(len(x[3]))
    temp = pandas.DataFrame(index = range(15),columns = [column + ' bes',column + ' ebs',column + ' ber',column + ' ebr'])
    temp[column+' bes'] = bes
    temp[column+' ebs'] = ebs
    temp[column+' ber'] = ber
    temp[column+' ebr'] = ebr
    total = pandas.concat([total,temp],axis = 1)

FVADat = data['FVA'].tolist()
changedRxns = list()
for x in FVADat:
    x = eval(x)
    types = list()
    for y in x:
        rxns = list()
        for z in y:
            rxns.append(z[0])
        types.append(rxns)
    changedRxns.append(types)
    print types
    print "###############\n\n"

uniqueChanged = list()
temp = pandas.DataFrame(index=range(15),columns=['UniqueFVAReact bes','UniqueFVAReact ebs','UniqueFVAReact ber','UniqueFVAReact ebr'])

for x,y in zip(changedRxns,range(15)):
    bes = list(set(x[0]) - set(x[1]))
    ebs = list(set(x[1]) - set(x[0]))
    ber = list(set(x[2]) - set(x[3]))
    ebr = list(set(x[3]) - set(x[2]))
    uniqueChanged.append([bes,ebs,ber,ebr])
    temp.loc[y,:] = [len(z) for z in uniqueChanged[-1]]

total = pandas.concat([total, temp], axis=1)

total.to_csv('NumberOfDifferences.csv')
#print total


#print KOReactions

#print data