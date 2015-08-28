# ---
# title: "make_residualized_case_dataset.R"
# author: "Kristen Dang"
# date: "Aug. 19, 2015"
# output: html_document
# ---

library('rGithubClient')
source('/Users/kristen/Computing/external_software/rgithubclient_authenticate.R')
sourceRepoFile(sageCode, 'scri_cran/cranio_common_routines.R')
setwd('~/Computing/cranio/')

cutoff = 0.05

## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)


## Get PALX20 protein-coding dataset with outliers removed
## Make case-only dataset
dataEnt = synGet('syn4738853')
load(getFileLocation(dataEnt))

metadataFiltered$Sample_Type = as.character(metadataFiltered$Sample_Type)
metadataFiltered$Sample_Type[grep("Coronal",metadataFiltered$Sample_Type)] = "Coronal"
metadataFiltered$Sample_Type[grep("Metopic",metadataFiltered$Sample_Type)] = "Metopic"
metadataFiltered$Sample_Type[grep("Sagittal",metadataFiltered$Sample_Type)] = "Sagittal"
metadataFiltered$Sample_Type = as.factor(metadataFiltered$Sample_Type)
minimalSetNoSex = c("Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES", "Initial_growth_duration_days")

caseMeta = metadataFiltered[-which(metadataFiltered$Sample_Type == "Control"),]
dim(caseMeta)
tempMeta = caseMeta[caseMeta$Px_Code %in% colnames(in.dge),]
dim(tempMeta)

## Run model and get residuals
tempMeta$Sample_Type = as.factor(as.character(tempMeta$Sample_Type))
tempModel = model.matrix(as.formula(paste("~",paste(minimalSetNoSex,collapse = "+"),sep = "")), data = tempMeta)

data.voom = voom(in.dge,tempModel,plot=FALSE) 
fit = lmFit(data.voom,tempModel)
save(fit,file = "resid_cases_pIGDD_fit.Robj.bz2",compress = "bzip2")
resid = residuals(fit,y = data.voom)
write.table(formatC(resid,digits=6,format="fg"),file = "scri-cran_resid_cases_pIGDD.tsv",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)

## Save data to Synapse
synStore(File(path = "scri-cran_resid_cases_pIGDD.tsv",parentId='syn4893931'))
synStore(File(path = "resid_cases_pIGDD_fit.Robj.bz2",parentId='syn4893931'))
