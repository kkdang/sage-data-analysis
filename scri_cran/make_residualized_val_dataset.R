# ---
# title: "make_residualized_dataset.R"
# author: "Kristen Dang"
# date: "Nov. 1, 2016"
# output: html_document
# ---

library(githubr)
library(edgeR)
source('/Users/kkdang/Computing/rgithubclient_authenticate.R')
sourceRepoFile(sageCode, "rnaseq_analysis_functions.R")
sourceRepoFile(sageCode, 'biomart_fxns.R')
setwd('~/Computing/cranio/')

cutoff = 0.05

## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata_validation.R")
metadataFiltered = processMetadataVal(plot = FALSE)
proteinCoding_genes = getByBiotype(biotype = 'protein_coding',gene = TRUE)


## Get data
dataEnt = read.delim(getFileLocation(synGet('syn7477102')), row.names = 1, header = TRUE)
lookupTable = synTableQuery('SELECT * FROM syn7477113')
head(lookupTable@values)

newNames = lookupTable@values$Investigator_Sample_Name[match(colnames(dataEnt), paste("X",lookupTable@values$Sample_Name, sep = ""))]
colnames(dataEnt) = newNames

data.dge = DGEList(counts=dataEnt,remove.zeros=TRUE)
palx20.dge = filterByFractionPresent(data.dge,fraction = 0.20)

## Remove samples that don't have seq data.
## Order the data in the same way for counts and metadata
toRemove = setdiff(rownames(metadataFiltered), colnames(palx20.dge))
metadataMatching = metadataFiltered[-which(rownames(metadataFiltered) %in% toRemove),-2]
metadataMatching = metadataMatching[match(colnames(palx20.dge), rownames(metadataMatching)),]
head(rownames(metadataMatching))
head(colnames(palx20.dge))
tail(rownames(metadataMatching))
tail(colnames(palx20.dge))

## Protein-coding genes, PALX=20%
pc.dge = palx20.dge[which(rownames(palx20.dge) %in% proteinCoding_genes[,1]),]
pc_palx20.dge = calcNormFactors(pc.dge)


## Run model and get residuals
tempModel = model.matrix(~Plate, data = metadataMatching)
data.voom = voom(pc_palx20.dge,tempModel,plot=FALSE) 
fit = lmFit(data.voom,tempModel)
#save(fit,file = "resid_pIGDD_fit.Robj.bz2",compress = "bzip2")
residVal = residuals(fit,y = data.voom)
#write.table(formatC(resid,digits=6,format="fg"),file = "resid_pIGDD.tsv",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)


## Save data to Synapse
#synStore(File(path = "resid_pIGDD.tsv",parentId='syn4893931'))
#synStore(File(path = "resid_pIGDD_fit.Robj.bz2",parentId='syn4893931'))


