from EssentialGeneFinder import *
from pandas import DataFrame
import cobra
import cobra.manipulation

def creationOfDeletedReaction():
    assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
    reference = {'GSE90620': [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1], 'GDS3572': [1, 1, 1, 0, 0, 0],
                 'GSE65870': [1, 1, 1, 1, 1, 1], 'GDS4244': [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                 'GSE30021': [0, 0, 0, 0, 0, 0, 1, 1, 1]}
    referenceMed = {'GSE90620': 1, 'GDS3572': 0, 'GSE65870': 2, 'GDS4244': 1, 'GSE30021': 1}
    nsamples = [12,12,6,9,6]
    deletions = [[set(),set()],[set(),set(),set()]]
    model = cobra.io.read_sbml_model("iPAE1146.xml")
    for i in range(len(assent)):
        res = findDeletedReactions(assent[i],model,nsamples[i],reference[assent[i]],referenceMed[assent[i]],'gurobi','C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')
        if not len(res[0]) == 0 : deletions[0][referenceMed[assent[i]]] = deletions[0][referenceMed[assent[i]]].union(res[0])
        else: print "None found in control: %d" % i
        if not len(res[1]) == 0 : deletions[1][referenceMed[assent[i]]] = deletions[1][referenceMed[assent[i]]].union(res[1])
        else: print "None found in exp: %d" % i

    return [deletions[0][1] & deletions[0][1],deletions[1][0] & deletions[1][1] & deletions[1][2]]
def findDeletedReactions(assent,CSmodel,numSamples,control,media,solver = 'gurobi',path= None):
    reactKOs = [set(), set()]

    for x in [z + 1 for z in range(numSamples)]:
        tempModel = CSmodel.copy()
        fileName = join(path, assent + '_' + str(x) + '.txt')
        statesFile = open(fileName, 'r')
        geneStates = statesFile.readlines()
        geneStates = [y.split(" ") for y in geneStates]
        geneNames = []
        for y in geneStates:
            if int(y[1]) == 0:
                geneNames.append(y[0])
        tempModel.solver = solver  # load model
        if not (len(geneNames) == 0):
            cobra.manipulation.delete_model_genes(tempModel, geneNames)
            temp = set()
            for r in tempModel.reactions:
                if r.bounds == (0,0):
                    temp.add(r.id)

            if len(reactKOs[control[x-1]]) == 0:
               reactKOs[control[x-1]] = temp.copy()
            else:
               reactKOs[control[x-1]] |= temp.copy()

    return reactKOs

data = creationOfDeletedReaction()
print data
file = open('removedReactionControl.txt','w')
file.write('Control\n')
[file.write(x+'\n') for x in data[0]]
file.close()
file = open('removedReactionExp.txt','w')
file.write('Experimental\n')
[file.write(x+'\n') for x in data[1]]
file.close()
