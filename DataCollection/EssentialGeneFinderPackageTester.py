from EssentialGeneFinder import *


def creationOfEssentialGeneData():
    assent = ['GSE90620','GDS4244','GDS3572','GSE30021','GSE65870']
    nsamples = [12,12,6,9,6]

    for i in range(5):

        createEssentialGeneDataAssent(assent[i],'.mat',nsamples[i],'C:\Users\Ethan Stancliffe\Desktop\Summer2017\Papin Lab\Pseudomonas Aeruginosa ABR Project\GeneEssentialityPA01')
      #  model = cobra.io.read_sbml_model("iPAE1146.xml")
       # createEssentialGeneModel(model,'IPAE1146')


#creationOfEssentialGeneData()

normal = ComparisionGene('IPAE1146','This is a test')


test = CSGene('GSE90620',[1,1,1,0,0,0,0,0,0,0,0,0],'LB')
test.processSamples()
#print test.sampleEG
#print test.findUniqueFrom(normal.getData())
print test.findUniqueFromControl()

