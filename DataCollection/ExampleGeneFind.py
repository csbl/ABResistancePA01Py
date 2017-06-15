import cobra.io
from EssentialGeneFinder import *


assent = 'GSE90620_'
nsamples = 12

createEssentialGeneDataAssent(assent,'.mat',nsamples,'C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\MATLAB_TIGER_Scripts')

#model = cobra.io.read_sbml_model("iPAE1146.xml")

#createEssentialGeneModel(model,'IPAE1146')
