# ---
# title: "make_residualized_dataset.R"
# author: "Kristen Dang"
# date: "Apr. 4, 2016"
# output: html_document
# ---

library(githubr)
library(edgeR)
source('/Users/kkdang/Computing/rgithubclient_authenticate.R')
sourceRepoFile(sageCode, "rnaseq_analysis_functions.R")
setwd('~/Computing/cranio/')

cutoff = 0.05

## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)
metadataFiltered$Sample_Type = as.character(metadataFiltered$Sample_Type)
metadataFiltered$Sample_Type[grep("Coronal",metadataFiltered$Sample_Type)] = "Coronal"
metadataFiltered$Sample_Type[grep("Metopic",metadataFiltered$Sample_Type)] = "Metopic"
metadataFiltered$Sample_Type[grep("Sagittal",metadataFiltered$Sample_Type)] = "Sagittal"
metadataFiltered$Sample_Type = as.factor(metadataFiltered$Sample_Type)

## Get counts data.
dataEnt = synGet('syn2820309')
estReads = read.csv(getFileLocation(dataEnt), row.names = 1)
colnames(estReads) = metadataFiltered$Px_Code[match(colnames(estReads), paste("X", metadataFiltered$SAMPLE,sep = ""))]
data.dge = DGEList(counts=estReads,remove.zeros=TRUE)
rm(estReads)

# Remove outliers and samples that don't have seq data.
outlierTable = synTableQuery('SELECT * FROM syn3354503')
outlierData = outlierTable@values
noOutliers.dge = data.dge[,-which(colnames(data.dge) %in% outlierData$minimalSet)]


## Order the data in the same way for counts and metadata
metadataMatching = metadataFiltered[-which(!metadataFiltered$Px_Code %in% colnames(noOutliers.dge)),]
metadataMatching = metadataMatching[match(colnames(noOutliers.dge), metadataMatching$Px_Code),]
head(metadataMatching$Px_Code)
head(colnames(noOutliers.dge))
dim(noOutliers.dge)

## Make PALX20 dataset 
palx20.dge = filterByFractionPresent(noOutliers.dge,fraction=0.1,minCount=2)
dim(palx20.dge)
palx20.dge = calcNormFactors(palx20.dge)


## Run model and get residuals
minimalSetNoSex = c("Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES", "Initial_growth_duration_days")
tempModel = model.matrix(as.formula(paste("~",paste(minimalSetNoSex,collapse = "+"),sep = "")), data = metadataMatching)
data.voom = voom(palx20.dge,tempModel,plot=FALSE) 
fit = lmFit(data.voom,tempModel)
#save(fit,file = "resid_pIGDD_fit.Robj.bz2",compress = "bzip2")
resid = residuals(fit,y = data.voom)
#write.table(formatC(resid,digits=6,format="fg"),file = "resid_pIGDD.tsv",quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)


## Save data to Synapse
#synStore(File(path = "resid_pIGDD.tsv",parentId='syn4893931'))
#synStore(File(path = "resid_pIGDD_fit.Robj.bz2",parentId='syn4893931'))


