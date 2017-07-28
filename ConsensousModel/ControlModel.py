
import cobra.manipulation

import DataCollection.EssentialGeneFinder

type = 'Control'
model = cobra.io.read_sbml_model("iPAE1146.xml")
DataCollection.EssentialGeneFinder.consensousModelCreator(model,type)

#cobra.manipulation.delete_model_genes(model,['PA4628'])
#res = cobra.flux_analysis.single_gene_deletion(model)
#res.to_csv(type+'GeneDDWPA4628')
