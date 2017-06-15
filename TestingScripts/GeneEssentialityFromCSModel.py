import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import cobra.io
import pandas
import numpy
from numpy import *


from bisect import bisect_left

def binary_search(a, x, lo=0, hi=None):  # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)
    pos = bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1)  # don't walk off the end  #https://stackoverflow.com/questions/212358/binary-search-bisection-in-python




CSmodel = cobra.io.load_matlab_model('CSCobraModel.mat')

CSmodel.solver = 'gurobi'

solution = CSmodel.optimize()

#print CSmodel.summary()

res = cobra.flux_analysis.single_gene_deletion(CSmodel)

#print res

essential = res['flux'] == 0

essGenes = res[essential]

#print essGenes

#print "Now to load the non context specific model\n"

model = model = cobra.io.read_sbml_model("iPAE1146.xml")
model.solver = 'gurobi'
sol = model.optimize()

res2 = cobra.flux_analysis.single_gene_deletion(model)


essential = res2['flux'] == 0

essGenes2 = res2[essential]

#print essGenes

genesNormal = sort(array(essGenes2.index.tolist()))
genesCS = sort(array(essGenes.index.tolist()))

#print genesNormal
#print genesCS

for x in genesCS:
    if binary_search(genesNormal,x) == -1:
        print x + 'is unique to the context specific model'
for x in genesNormal:
    if binary_search(genesCS,x) == -1:
        print x + 'is unique to the general model'



