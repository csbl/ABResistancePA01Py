import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import cobra.io
from numpy import *
from os.path import join

def createEssentialGeneModelDD(model,label,type = 'smbl',solver = 'gurobi'):
    model.solver = solver
    res = cobra.flux_analysis.double_gene_deletion(model)
    newFile = open(('EssGene'+label+'.txt'),'w')
    res.sort_index()
    print res


model = cobra.io.read_sbml_model('IPAE1146.xml')

model.optimize()

createEssentialGeneModelDD(model,'IPAE1146')

