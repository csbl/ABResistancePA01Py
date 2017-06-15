

import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import matplotlib
import numpy
import matlab.engine
import matlab
import scipy.io

geneExpression = scipy.io.loadmat('GeneExpresssion.mat')

print(geneExpression)
model = cobra.io.read_sbml_model("iPAE1146.xml")

print model.genes

solution = model.optimize() #optimize
