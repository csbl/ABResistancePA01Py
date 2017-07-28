import DataCollection.EssentialGeneFinder
import cobra.io
import ModelComparison
from multiprocessing import Process
from pandas import DataFrame,concat

if __name__ == '__main__':
    model = cobra.io.read_sbml_model('iPAE1146.xml')
    geneListF = open('sharedGenes.txt', 'r')
    geneListShared = [x[:-1] for x in geneListF.readlines()]
    geneListF.close()
    geneListF = open('RUniqGenes.txt', 'r')
    geneListR = [x[:-1] for x in geneListF.readlines()]
    geneListF.close()
    geneListF = open('CUniqGenes.txt', 'r')
    geneListC = [x[:-1] for x in geneListF.readlines()]
    geneListF.close()


    results = dict()
    type = {1: 'R', 0: 'C'}
    orderOfDel = [[geneListC, []]] + [[geneListC, [x]] for x in geneListR] + [[[], [x]] for x in geneListR] + [[[x], []] for x in geneListC]
    mapDel2Res = {x:y for x,y in zip(range(len(orderOfDel)),orderOfDel)}

    i = 0
    for x in orderOfDel:

        p = Process(target = DataCollection.EssentialGeneFinder.consensousModelCreator, args = (model.copy(),'Resistant',list(geneListShared) + list(x[1])))
        q = Process(target = DataCollection.EssentialGeneFinder.consensousModelCreator,args=(model.copy(),'Control',list(geneListShared)+list(x[0])))
        p.start()
        q.start()
        p.join()
        q.join()
        results[i] = ModelComparison.run()
        i += 1

    data = DataFrame.from_dict(results,orient='index')
    data = concat([data,DataFrame.from_dict(mapDel2Res,orient = 'index')],axis = 1)
    print data




