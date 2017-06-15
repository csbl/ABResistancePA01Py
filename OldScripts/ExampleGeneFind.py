import cobra.io
from DataCollection.EssentialGeneFinder import *


assent = 'GSE90620_'
nsamples = 12

createEssentialGeneDataAssent(assent,'.mat',nsamples,'C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')

#model = cobra.io.read_sbml_model("iPAE1146.xml")

#createEssentialGeneModel(model,'IPAE1146')
