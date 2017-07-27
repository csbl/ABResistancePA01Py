import cobra
import cobra.flux_analysis
import cobra.solvers.gurobi_solver
import cobra.io
from numpy import *
from os.path import join
from scipy.stats import *
from bisect import bisect_left
import cobra.manipulation
from copy import deepcopy
import scipy.stats
import DataCollection.EssentialGeneFinder

type = 'Resistant'
model = cobra.io.read_sbml_model("iPAE1146.xml")
geneListF = open(type+'Genes.txt','r')
geneList = geneListF.readlines()
geneList = [x[:-1] for x in geneList]
cobra.manipulation.delete_model_genes(model,geneList)
removedReactionsID = []
for x in model.reactions:
    if x.bounds == (0,0):
        removedReactionsID.append(x.id)
geneListF.close()
outputFile = open(type +'ReactionsKO.txt','w')
#outputFile.write('ReactionID ReactionName')
[outputFile.write(x+'\n') for x in removedReactionsID]
outputFile.close()

res = cobra.flux_analysis.single_gene_deletion(model)
res.to_csv(type+'GeneDeletion.csv')
res = cobra.flux_analysis.single_reaction_deletion(model)
res.to_csv(type + 'RxnDeletion.csv')
res = cobra.flux_analysis.flux_variability_analysis(model)
res.to_csv(type + 'FVA.csv')
