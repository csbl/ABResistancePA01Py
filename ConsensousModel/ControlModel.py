import cobra.manipulation
import CSModelDataScripts.ABResistanceAnalysis
"""Script for creating control consensous model using the consensousModelCreator function in the ABResistance package"""
type = 'Control'
model = cobra.io.read_sbml_model("iPAE1146.xml")
CSModelDataScripts.ABResistanceAnalysis.consensusModelCreator(model, type)
