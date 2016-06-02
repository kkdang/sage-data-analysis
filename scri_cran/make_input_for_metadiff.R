# make input file for meta-diff
# Dec. 29, 2015
# KKD for Sage Bionetworks


library(synapseClient)
library('rGithubClient')
source('/Users/kristen/Computing/external_software/rgithubclient_authenticate.R')
setwd('~/Computing/cranio/')
synapseLogin()

## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)
metadataFiltered$Sample_Type = as.character(metadataFiltered$Sample_Type)
metadataFiltered$Sample_Type[grep("Coronal",metadataFiltered$Sample_Type)] = "Coronal"
metadataFiltered$Sample_Type[grep("Metopic",metadataFiltered$Sample_Type)] = "Metopic"
metadataFiltered$Sample_Type[grep("Sagittal",metadataFiltered$Sample_Type)] = "Sagittal"
metadataFiltered$Sample_Type = as.factor(metadataFiltered$Sample_Type)

## Remove outliers
outlierTable = synTableQuery('SELECT * FROM syn3354503')
outlierData = outlierTable@values
noOutliers.meta = metadataFiltered[-which(metadataFiltered$Px_Code %in% outlierData$transcriptsSet),]


## Restrict to model variables, set column names and order as required by meta-diff
modelToKeep = c("Age_mos.","PCT_CORRECT_STRAND_READS","Initial_growth_duration_days", "SAMPLE", "caseStatus", "Px_Code")
colnames(noOutliers.meta)
model.meta = noOutliers.meta[,which(colnames(noOutliers.meta) %in% modelToKeep)]
head(model.meta)
model.meta$"File_Name" = paste("/mnt/", model.meta$SAMPLE, ".isoforms.fpkm_tracking", sep = "")
model.meta$"C_group" = rep(0,nrow(model.meta))
model.meta$C_group[which(model.meta$caseStatus == "case")] = 1
model.meta$Sample = model.meta$Px_Code
tail(model.meta)

toRemove = c("caseStatus", "SAMPLE", "Px_Code")
final.model = model.meta[,-which(colnames(model.meta) %in% toRemove)]
head(final.model)
toOrder = c("Sample", "File_Name", "C_group", "Age_mos.", "PCT_CORRECT_STRAND_READS", "Initial_growth_duration_days")
output.file = final.model[,match(toOrder, colnames(final.model))]
head(output.file)
formatted.output = output.file
formatted.output$Age_mos. = formatC(output.file$Age_mos., digits=6, format="fg")
formatted.output$Initial_growth_duration_days = formatC(output.file$Initial_growth_duration_days, digits=6, format="fg")
formatted.output$Age_mos. = formatC(output.file$Age_mos., digits=6, format="fg")
write.table(formatted.output, file = "metaDiff_caseControl_input.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

toRemoveNA = grep(pattern = 'NA', x = formatted.output$File_Name)
formattedShort.output = formatted.output[-toRemoveNA,]


table(metadataFiltered$Sample_Type)
coronalSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Coronal"])
controlSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Control"])
lambdoidSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Lambdoid"])
metopicSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Metopic"])
sagittalSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Sagittal"])

## Write out file for specific suture-types vs control
toWrite = which(formattedShort.output$Sample %in% c(coronalSampleIds, controlSampleIds))
write.table(formattedShort.output[toWrite,], file = "metaDiff_coronalControl_input.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


toWrite = which(formattedShort.output$Sample %in% c(metopicSampleIds, controlSampleIds))
write.table(formattedShort.output[toWrite,], file = "metaDiff_metopicControl_input.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

toWrite = which(formattedShort.output$Sample %in% c(sagittalSampleIds, controlSampleIds))
write.table(formattedShort.output[toWrite,], file = "metaDiff_sagittalControl_input.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
