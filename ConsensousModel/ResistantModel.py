
import cobra.manipulation
import DataCollection.EssentialGeneFinder

type = 'Resistant'
model = cobra.io.read_sbml_model("iPAE1146.xml")
DataCollection.EssentialGeneFinder.consensousModelCreator(model,type)


#res = cobra.flux_analysis.single_gene_deletion(model)
#res.to_csv(type+'GeneDDWPA4628')
