

import cobra
import os
from os.path import join
import cobra.flux_analysis
import cobra.solvers.gurobi_solver



#model = cobra.test.create_test_model("textbook")
#model = cobra.io.read_sbml_model(join(data_dir,"iPAE1146.xml")) #load model
model = cobra.io.read_sbml_model("iPAE1146.xml")
solution = model.optimize() #optimize
print(solution)

print model.summary()

sol = cobra.flux_analysis.single_gene_deletion(model)
sol = sol[sol['flux'] > 1]
print sol[sol['flux'] < 15.5]
#geneExpression = eng.testRun()

#tissueModel = eng.testTissueModel()

#cobra.io.save_matlab_model(model,"myModel")

#tissueModel = eng.createTissueSpecificModel(eng.load('myModel.mat'),geneExpression)
#print(len(model.reactions))
#print(len(model.metabolites))
#print(len(model.genes))

#tissueModel.summary()

#print model.genes


#delResults = cobra.flux_analysis.single_gene_deletion(model)#Attmept to simulate a gene deletion.. maybe
#print(delResults)



#model.solver = 'gurobi'
#print(type(model.solver))
#print(cobra.solvers.solver_dict)
#solver = cobra.solvers.gurobi_solver
#lp = solver.create_problem(model,objective_sense="maximize") #Check that solver is working
#print(solver.solve_problem(lp))