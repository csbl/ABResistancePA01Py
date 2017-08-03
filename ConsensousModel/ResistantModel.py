
import cobra.manipulation
import CSModelDataScripts.ABResistanceAnalysis

"""Script for creating resistant consensous models using the consensousModelCreator function in the ABResistanceAnalysis package"""
type = 'Resistant'
model = cobra.io.read_sbml_model("iPAE1146.xml")
CSModelDataScripts.ABResistanceAnalysis.consensusModelCreator(model, type)


