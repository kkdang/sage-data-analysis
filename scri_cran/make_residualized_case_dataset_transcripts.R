# ---
# title: "make_residualized_case_dataset.R"
# author: "Kristen Dang"
# date: "Sept. 28, 2015"
# output: html_document
# ---

library('rGithubClient')
source('/Users/kristen/Computing/external_software/rgithubclient_authenticate.R')
sourceRepoFile(sageCode, 'scri_cran/cranio_common_routines.R')
setwd('~/Computing/cranio/')
proteinCoding_transcripts = getByBiotype(gene = FALSE)

## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)
metadataFiltered$Sample_Type = as.character(metadataFiltered$Sample_Type)
metadataFiltered$Sample_Type[grep("Coronal",metadataFiltered$Sample_Type)] = "Coronal"
metadataFiltered$Sample_Type[grep("Metopic",metadataFiltered$Sample_Type)] = "Metopic"
metadataFiltered$Sample_Type[grep("Sagittal",metadataFiltered$Sample_Type)] = "Sagittal"
metadataFiltered$Sample_Type = as.factor(metadataFiltered$Sample_Type)

## Remove outliers
data.dge = generateDataObj('syn3064954')
outlierTable = synTableQuery('SELECT * FROM syn3354503')
outlierData = outlierTable@values
noOutliers.dge = data.dge[,-which(colnames(data.dge) %in% outlierData$transcriptsSet)]

## Remove samples that don't have seq data.
## Order the data in the same way for counts and metadata
metadataMatching = metadataFiltered[-which(!metadataFiltered$Px_Code %in% colnames(noOutliers.dge)),]
metadataMatching = metadataMatching[match(colnames(noOutliers.dge), metadataMatching$Px_Code),]
head(metadataMatching$Px_Code)
head(colnames(noOutliers.dge))


## PALX=5%, cases only
palx5.dge = filterByFractionPresent(noOutliers.dge,fraction=0.05,minCount=3)
palx5.dge = calcNormFactors(palx5.dge)
cases.dge = palx5.dge[,-which(metadataMatching$Sample_Type == "Control")]
dim(cases.dge)
#save(cases.dge,file = "transcripts_cases_DGE.Robj.bz2",compress = "bzip2")
#synStore(File(path = "transcripts_cases_DGE.Robj.bz2",parentId='syn2820780'),forceVersion=FALSE)


minimalSetNoSex = c("Age_mos.","PCT_CORRECT_STRAND_READS","Initial_growth_duration_days")

caseMeta = metadataMatching[-which(metadataMatching$Sample_Type == "Control"),]
dim(caseMeta)
head(colnames(cases.dge))
head(caseMeta$Px_Code)
tail(colnames(cases.dge))
tail(caseMeta$Px_Code)

## Run model and get residuals
caseMeta$Sample_Type = as.factor(as.character(caseMeta$Sample_Type))
tempModel = model.matrix(as.formula(paste("~",paste(minimalSetNoSex,collapse = "+"),sep = "")), data = caseMeta)

data.voom = voom(cases.dge,tempModel,plot=FALSE) 
fit = lmFit(data.voom,tempModel)
#save(fit,file = "transcripts_resid_cases_pIGDD_fit.Robj.bz2",compress = "bzip2")
resid = residuals(fit,y = data.voom)
#write.table(formatC(resid,digits=6,format="fg"),file = "transcripts_resid_cases_pIGDD.tsv",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)

## Save data to Synapse
#synStore(File(path = "transcripts_resid_cases_pIGDD.tsv",parentId='syn4893931'))
#synStore(File(path = "transcripts_resid_cases_pIGDD_fit.Robj.bz2",parentId='syn4893931'))
