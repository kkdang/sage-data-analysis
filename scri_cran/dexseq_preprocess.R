#! /usr/bin/env Rscript
# title: "dexseq_preprocess.Rmd"
# author: "Kristen Dang"
# date: "June 20, 2016"

library(synapseClient)
synapseLogin()
library(githubr)
source('~/Music/rgithubclient_authenticate.R')
sourceRepoFile(sageCode, "rnaseq_analysis_functions.R")
library(DEXSeq)
source('/work/DAT_112__Craniosynostosis_seq_subtype/Scripts/load_SubreadOutput.R')
setwd('/work/DAT_112__Craniosynostosis_seq_subtype/Results/')


## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)
metadataFiltered$Sample_Type = as.character(metadataFiltered$Sample_Type)
metadataFiltered$Sample_Type[grep("Coronal",metadataFiltered$Sample_Type)] = "Coronal"
metadataFiltered$Sample_Type[grep("Metopic",metadataFiltered$Sample_Type)] = "Metopic"
metadataFiltered$Sample_Type[grep("Sagittal",metadataFiltered$Sample_Type)] = "Sagittal"
metadataFiltered$Sample_Type = as.factor(metadataFiltered$Sample_Type)


## Get exon data and read into dexseq object
exonCountsEnt = synGet('syn5568207')
exonCounts = read.delim(getFileLocation(exonCountsEnt), header=TRUE,stringsAsFactors = FALSE,row.names = 1)
gtfEnt = synGet('syn5386238')

colnames(exonCounts)
x = metadataFiltered$Px_Code[match(colnames(exonCounts), paste("X",metadataFiltered$SAMPLE, sep = ""))]
colnames(exonCounts) = x



metadataMatching = metadataFiltered[which(metadataFiltered$Px_Code %in% colnames(exonCounts)),]
rownames(metadataMatching) = metadataMatching$Px_Code
colnames(metadataMatching)
colnames(metadataMatching)[which(colnames(metadataMatching) == "caseStatus")] = "condition"
metadataFormatted = metadataMatching[,-1]
head(metadataFormatted)
levels(metadataFormatted$condition)
metadataFormatted$condition = relevel(metadataFormatted$condition,ref = "control")
levels(metadataFormatted$condition)


# Remove genes which are not detected in any sample.
# Remove one more set of genes that are not in the gtf.
allZero = which(rowSums(exonCounts) == 0)
zeroPurged = exonCounts[-allZero,]
toRemove = "ENSG00000261934+ENSG00000253910+ENSG00000204956+ENSG00000253159+ENSG00000242419+ENSG00000253846+ENSG00000262209+ENSG00000253485+ENSG00000253767+ENSG00000240184+ENSG00000240764+ENSG00000253731+ENSG00000254221+ENSG00000081853+ENSG00000253305+ENSG00000254122"
splitNames = sapply(as.list(rownames(zeroPurged)),function(x) {unlist(strsplit(x = x,split = ":"))[1]})
head(splitNames)
tail(splitNames)
x = which(splitNames %in% toRemove)
palo = zeroPurged[-x,]


splitNames = sapply(as.list(rownames(palo)),function(x) {unlist(strsplit(x = x,split = ":"))[1]})
x = which(splitNames %in% toRemove)


# PALX filtering
palx80 = filterByFractionPresent(palo,fraction = 0.8)

data.dex = DEXSeqDataSetFromFeatureCounts(countData = palx80,sampleData = metadataFormatted ,flattenedfile = getFileLocation(gtfEnt))
head(data.dex)
data.dex = estimateSizeFactors(data.dex)
save(data.dex,file = "data_wSizeFactors.Robj.bz",compress = "bzip2")

rm(exonCounts)
gc()


## QC analysis
# Density plots
# allCounts.norm = log2(featureCounts(data.dex,normalized=TRUE))
# dim(allCounts.norm)
# head(allCounts.norm)
# plot(density(allCounts.norm[,1]), col = "blue", main = "norm data", ylim = range(0,1))
# for (i in 2:ncol(allCounts.norm)) { lines(density(allCounts.norm[,i]), col = "blue") }
# 
# # Median - IQR
# meds = apply(allCounts.norm,MARGIN = 2,FUN = median,na.rm = TRUE)
# iqrs = apply(allCounts.norm,MARGIN = 2,FUN = IQR, na.rm = TRUE)
# plot(meds,iqrs, xlab = "medians", ylab = "interquartile ranges", main = "log2(CPM) protein-coding genes")
# abline(h = 8.6, col = "sandybrown", lty = 2)
# abline(v = 2.98, col = "sandybrown", lty = 2)
# which(iqrs > 8.6)
# which(meds < 2.98)


