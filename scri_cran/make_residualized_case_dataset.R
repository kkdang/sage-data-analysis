# ---
# title: "make_residualized_case_dataset.R"
# author: "Kristen Dang"
# date: "June 18, 2015"
# output: html_document
# ---
  
library('rGithubClient')
token = read.delim('~/Movies/rGHclient_token.txt',header = FALSE)
setGithubToken(as.character(token[1,1]))
sageCode = getRepo(repository="kkdang/sage-data-analysis")
sourceRepoFile(sageCode, 'scri_cran/cranio_common_routines.R')
proteinCoding_genes = getByBiotype()
setwd('~/Computing/cranio/')

cutoff = 0.05

## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)


## Get PALX20 protein-coding dataset with outliers removed
data.dge = generateDataObj('syn2820309')

outlierTable = synTableQuery('SELECT * FROM syn3354503')
outlierData = outlierTable@values
noOutliers.dge = data.dge[,-which(colnames(data.dge) %in% outlierData$minimalSet)]

## Remove samples that don't have seq data.
## Order the data in the same way for counts and metadata
metadataMatching = metadataFiltered[-which(!metadataFiltered$Px_Code %in% colnames(noOutliers.dge)),]
metadataMatching = metadataMatching[match(colnames(noOutliers.dge), metadataMatching$Px_Code),]
head(metadataMatching$Px_Code)
head(colnames(noOutliers.dge))

## Protein-coding genes, PALX=20%
b = as.list(rownames(noOutliers.dge))
strippedNames = sapply(b,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})  
pc.dge = noOutliers.dge[which(strippedNames %in% proteinCoding_genes[,1]),]
pc_palo.dge = DGEList(counts = getCounts(pc.dge),group = pc.dge$samples$group)
pc_palx20.dge = filterByFractionPresent(pc_palo.dge,fraction=0.20,minCount=3)
pc_palx20.dge = calcNormFactors(pc_palx20.dge)

## Make case-only dataset
in.dge = pc_palx20.dge[,-which(metadataMatching$Sample_Type == "Control")]
dim(in.dge)
save(in.dge,file = "pc_cases_DGE.Robj.bz2",compress = "bzip2")
synStore(File(path = "pc_cases_DGE.Robj.bz2",parentId='syn2820780'))


metadataMatching$Sample_Type = as.character(metadataMatching$Sample_Type)
metadataMatching$Sample_Type[grep("Coronal",metadataMatching$Sample_Type)] = "Coronal"
metadataMatching$Sample_Type[grep("Metopic",metadataMatching$Sample_Type)] = "Metopic"
metadataMatching$Sample_Type[grep("Sagittal",metadataMatching$Sample_Type)] = "Sagittal"
metadataMatching$Sample_Type = as.factor(metadataMatching$Sample_Type)
minimalSet = c("Sex","Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES")

caseMeta = metadataMatching[-which(metadataMatching$Sample_Type == "Control"),]
dim(caseMeta)

## Run model and get residuals
### Using residuals of model that does NOT include AlkP or sample_type
plotCols = c("grey63","white", "coral","gold2")
tempMeta = caseMeta
tempMeta$Sample_Type = as.factor(as.character(tempMeta$Sample_Type))
tempModel = model.matrix(as.formula(paste("~",paste(minimalSet,collapse = "+"),sep = "")), data = tempMeta)
data.voom = voom(in.dge,tempModel,plot=FALSE) 
fit = lmFit(data.voom,tempModel)
save(fit,file = "resid_cases_fit.Robj.bz2",compress = "bzip2")
resid = residuals(fit,y = data.voom)
write.table(formatC(resid,digits=6,format="fg"),file = "scri-cran_resid_cases.tsv",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)

## Save data to Synapse
synStore(File(path = "scri-cran_resid_cases.tsv",parentId='syn2820780'))
synStore(File(path = "resid_cases_fit.Robj.bz2",parentId='syn2820780'))
