import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import cobra.io
from numpy import *
from os.path import join
from pandas import DataFrame

def grabIndexValue(coordinates,dataFrame):
    x = dataFrame.index[coordinates[0]]
    y = dataFrame.columns.values[coordinates[1]]
    index = str(x) + ' ' +  str(y)
    value = array(dataFrame)[coordinates[0],coordinates[1]]
    return index, value


def createEssentialGeneModelDD(model,label,type = 'smbl',solver = 'gurobi'):
    model.solver = solver
    #res = cobra.flux_analysis.double_gene_deletion(model,return_frame = True ,number_of_processes=1)
    #file = open('DDTest.csv','w')
    newFile = open(('EssGene' + label + '.txt'), 'w')
    res = DataFrame.from_csv('DDTest.csv')
    xLabel = res.index
    yLabel = res.columns.values
    res = array(res)
    for i in range(len(xLabel)):
        for j in range(len(yLabel)):
            x = xLabel[i]
            y = yLabel[j]
            index = str(x) + '_' + str(y)
            value = res[i,j]
            newFile.write(index + ' ' + str(value) + '\n')
    newFile.close()
    #file.close


model = cobra.io.read_sbml_model('IPAE1146.xml')

model.optimize()

createEssentialGeneModelDD(model,'IPAE1146DD')

