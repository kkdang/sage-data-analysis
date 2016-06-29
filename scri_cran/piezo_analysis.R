#! /usr/bin/env Rscript
# KKD for Sage Bionetworks
# Jun 27, 2016
# Plot expression of PIEZO1

library(synapseClient)
synapseLogin()
library("VennDiagram")
library('gplots')
library('RColorBrewer')
library(ReactomePA)
library('githubr')
source('/Users/kkdang/Computing/rgithubclient_authenticate.R')
sourceRepoFile(sageCode, 'scri_cran/cranio_common_routines.R')
sourceRepoFile(sageCode, 'biomart_fxns.R')
setwd('~/Computing/cranio/')

## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)
metadataFiltered$Sample_Type = as.character(metadataFiltered$Sample_Type)
metadataFiltered$Sample_Type[grep("Coronal",metadataFiltered$Sample_Type)] = "Coronal"
metadataFiltered$Sample_Type[grep("Metopic",metadataFiltered$Sample_Type)] = "Metopic"
metadataFiltered$Sample_Type[grep("Sagittal",metadataFiltered$Sample_Type)] = "Sagittal"
metadataFiltered$Sample_Type = as.factor(metadataFiltered$Sample_Type)
head(metadataFiltered)


# Get groups of interest
groups = read.delim("~/Computing/cranio/piezo1_sample_groups.txt", header = TRUE)
groups
controls = which(metadataFiltered$caseStatus == "control")
cSamps = as.data.frame(cbind(metadataFiltered$SAMPLE[controls], rep(x = "control",times=length(controls))))
colnames(cSamps) = c("SAMPLE", "PIEZOgroup")
allSamples = rbind(groups, cSamps)


#Get dataset, protein-coding genes, PALX=20%, with minimal model
sourceRepoFile(sageCode, 'scri_cran/make_residualized_dataset.R')
proteinCoding_genes = getByBiotype()
cutoff = 0.3
b = as.list(rownames(noOutliers.dge))
strippedNames = sapply(b,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})  
pc.dge = noOutliers.dge[which(strippedNames %in% proteinCoding_genes[,1]),]
pc_palo.dge = DGEList(counts = getCounts(pc.dge),group = pc.dge$samples$group)
pc_palx50.dge = filterByFractionPresent(pc_palo.dge,fraction=0.50,minCount=3)
pc_palx50.dge = calcNormFactors(pc_palx50.dge)
rm(pc.dge)


#Functions
computeFit=function(in.dge,design,inContrast=my.contrasts[,1],plotTitle=modelName,pval=cutoff){
  data.voom = voom(in.dge,design,plot=FALSE)
  fit = lmFit(data.voom,design)
  testResults = calculateDE(in.fit=fit,inContrast=inContrast,pcutoff=pval,plotTitle=plotTitle,limma=TRUE)
  DEgeneTable = topTable(testResults$fit,number = nrow(in.dge))
  xSig = rownames(DEgeneTable)[DEgeneTable$adj.P.Val < pval]
  xStatus = as.numeric(rownames(testResults$fit) %in% xSig)
  limma::plotMA(testResults$fit,status = xStatus,cex=c(0.5),legend = FALSE,col=c("green3"),main = modelName)
  print(length(which(DEgeneTable$adj.P.Val<pval)))
  return(DEgeneTable)
}

# Filter metadata, add Piezo1 groups
shortMeta = metadataFiltered[-which(is.na(metadataFiltered$SAMPLE)),]
temp = merge(shortMeta, allSamples, by.x = "SAMPLE", by.y = "SAMPLE", all.x = FALSE, all.y = FALSE)
PiezoMeta = temp[temp$Px_Code %in% colnames(pc_palx50.dge),]
head(PiezoMeta)
duplicated(PiezoMeta$SAMPLE)
duplicated(PiezoMeta$Px_Code)
dim(PiezoMeta)
rm(temp)

# All groups v control signatures 
piezoData.dge = pc_palx50.dge[,colnames(pc_palx50.dge) %in% PiezoMeta$Px_Code]
dim(piezoData.dge)
PiezoMeta = PiezoMeta[match(PiezoMeta$Px_Code, colnames(piezoData.dge)),]
# Look for any genes without expression in these samples or without expression in the cases
which(rowSums(getCounts(piezoData.dge)) == 0)
cases = which(PiezoMeta$caseStatus == "case")
which(rowSums(getCounts(piezoData.dge)[,cases]) == 0)
x = rowSums(getCounts(piezoData.dge)[,cases]) 

# Remove genes with low average expression
hist(cpm(piezoData.dge,normalized.lib.sizes = TRUE,log = TRUE))
piezoAve = aveLogCPM(piezoData.dge,normalized.lib.sizes = TRUE)
hist(piezoAve)
piezoDataFiltered.dge = piezoData.dge[-which(piezoAve < 0),]
dim(piezoDataFiltered.dge)

head(colnames(piezoDataFiltered.dge))
head(PiezoMeta$Px_Code)

modelName = "case-control + full"
minimalSet = c("Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES", "Initial_growth_duration_days")
ccModel = model.matrix(as.formula(paste("~0+caseStatus",paste(minimalSet,collapse = "+"),sep = "+")), data = PiezoMeta)
rownames(ccModel) = colnames(piezoDataFiltered.dge)
my.contrasts = makeContrasts(caseVScontrol = caseStatuscase-caseStatuscontrol, levels=ccModel)  
results_ccFull = computeFit(in.dge = piezoDataFiltered.dge,design = ccModel)
de_all = rownames(results_ccFull)[which(results_ccFull$adj.P.Val<cutoff)]



# Groups 1-2 v control
dim(PiezoMeta)
dim(piezoData.dge)
groups12 = PiezoMeta$Px_Code[PiezoMeta$PIEZOgroup %in% c("1a", "1b", "2", "control")]
piezoData.dge = piezoData.dge[,colnames(piezoData.dge) %in% groups12]
dim(piezoData.dge)
g12PiezoMeta = PiezoMeta[PiezoMeta$Px_Code %in% colnames(piezoData.dge),]
g12PiezoMeta = g12PiezoMeta[na.omit(match(g12PiezoMeta$Px_Code, colnames(piezoData.dge))),]
head(g12PiezoMeta$Px_Code)
head(colnames(piezoData.dge))
# Look for any genes without expression in these samples or without expression in the cases
which(rowSums(getCounts(piezoData.dge)) == 0)
cases = which(g12PiezoMeta$caseStatus == "case")
which(rowSums(getCounts(piezoData.dge)[,cases]) == 0)

# Remove genes with low average expression
piezoDataFiltered.dge = piezoData.dge[-which(piezoAve < 0),]
dim(piezoDataFiltered.dge)


modelName = "case-control Piezo groups 1 & 2"
minimalSet = c("Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES", "Initial_growth_duration_days")
ccModel = model.matrix(as.formula(paste("~0+caseStatus",paste(minimalSet,collapse = "+"),sep = "+")), data = g12PiezoMeta)
rownames(ccModel) = colnames(piezoDataFiltered.dge)
my.contrasts = makeContrasts(caseVScontrol = caseStatuscase-caseStatuscontrol, levels=ccModel)  
results_ccG12 = computeFit(in.dge = piezoDataFiltered.dge,design = ccModel)
de_allG12 = rownames(results_ccG12)[which(results_ccG12$adj.P.Val<cutoff)]



# Group 1 v control
dim(PiezoMeta)
dim(pc_palx50.dge)
piezoData.dge = pc_palx50.dge[,colnames(pc_palx50.dge) %in% PiezoMeta$Px_Code]
dim(piezoData.dge)
PiezoMeta = PiezoMeta[match(PiezoMeta$Px_Code, colnames(piezoData.dge)),]

group1 = PiezoMeta$Px_Code[PiezoMeta$PIEZOgroup %in% c("1a", "1b", "control")]
piezoDataG1.dge = piezoData.dge[,colnames(piezoData.dge) %in% group1]
dim(piezoDataG1.dge)
g1PiezoMeta = PiezoMeta[PiezoMeta$Px_Code %in% colnames(piezoDataG1.dge),]
g1PiezoMeta = g1PiezoMeta[na.omit(match(g1PiezoMeta$Px_Code, colnames(piezoDataG1.dge))),]
head(g1PiezoMeta$Px_Code)
head(colnames(piezoDataG1.dge))
# Look for any genes without expression in these samples or without expression in the cases
which(rowSums(getCounts(piezoDataG1.dge)) == 0)
cases = which(g1PiezoMeta$caseStatus == "case")
which(rowSums(getCounts(piezoDataG1.dge)[,cases]) == 0)

# Remove genes with low average expression
piezoDataFiltered.dge = piezoDataG1.dge[-which(piezoAve < 0),]
dim(piezoDataFiltered.dge)


modelName = "case-control Piezo groups 1"
minimalSet = c("Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES", "Initial_growth_duration_days")
ccModel = model.matrix(as.formula(paste("~0+caseStatus",paste(minimalSet,collapse = "+"),sep = "+")), data = g1PiezoMeta)
rownames(ccModel) = colnames(piezoDataFiltered.dge)
my.contrasts = makeContrasts(caseVScontrol = caseStatuscase-caseStatuscontrol, levels=ccModel)  
results_ccG1 = computeFit(in.dge = piezoDataFiltered.dge,design = ccModel)




# Group 1a v 1b
dim(PiezoMeta)
dim(piezoData.dge)

group1 = PiezoMeta$Px_Code[PiezoMeta$PIEZOgroup %in% c("1a", "1b")]
piezoDataG1.dge = piezoData.dge[,colnames(piezoData.dge) %in% group1]
dim(piezoDataG1.dge)
g1PiezoMeta = PiezoMeta[PiezoMeta$Px_Code %in% colnames(piezoDataG1.dge),]
g1PiezoMeta = g1PiezoMeta[na.omit(match(g1PiezoMeta$Px_Code, colnames(piezoDataG1.dge))),]
g1PiezoMeta$PIEZOgroup = as.factor(g1PiezoMeta$PIEZOgroup)
# Look for any genes without expression in these samples or without expression in the cases
which(rowSums(getCounts(piezoDataG1.dge)) == 0)
cases = which(g1PiezoMeta$caseStatus == "case")
which(rowSums(getCounts(piezoDataG1.dge)[,cases]) == 0)

# Remove genes with low average expression
piezoDataFiltered.dge = piezoDataG1.dge[-which(piezoAve < 0),]
dim(piezoDataFiltered.dge)


modelName = "case-control Piezo groups 1a v 1b"
minimalSet = c("Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES", "Initial_growth_duration_days")
ccModel = model.matrix(as.formula(paste("~0+PIEZOgroup",paste(minimalSet,collapse = "+"),sep = "+")), data = g1PiezoMeta)
rownames(ccModel) = colnames(piezoDataFiltered.dge)
my.contrasts = makeContrasts(aVSb = PIEZOgroup1b-PIEZOgroup1a, levels=ccModel)  
results_ccG1ab = computeFit(in.dge = piezoDataFiltered.dge,design = ccModel)








# Using a model that also controls for sex
sourceRepoFile(sageCode, 'scri_cran/make_residualized_dataset.R')

msSex = c("Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES", "Initial_growth_duration_days", "Sex")
tempModel = model.matrix(as.formula(paste("~",paste(msSex,collapse = "+"),sep = "")), data = metadataMatching)
data.voom = voom(palx20.dge,tempModel,plot=FALSE) 
fit = lmFit(data.voom,tempModel)
altData = residuals(fit,y = data.voom)

piezoData  = cbind(altData[which(rownames(altData ) == "ENSG00000103335"),], gsub(pattern = "X",replacement = "",x = colnames(altData)))
colnames(piezoData) = c("ENSG00000103335", "Px_Code")
piezoData = as.data.frame(piezoData)
piezoData$ENSG00000103335 = as.numeric(as.character(piezoData$ENSG00000103335))
head(piezoData)
head(PiezoMeta)
PiezoMerged = merge(PiezoMeta, piezoData, by.x = "Px_Code", by.y = "Px_Code", all.x = TRUE, all.y = FALSE)
head(PiezoMerged)
beeswarm(PiezoMerged$ENSG00000103335~PiezoMerged$PIEZOgroup, ylab = "Piezo1 residualized exp")



# Original model, but including outlier sample 2088 (aka 95689)
outlierData = outlierTable@values
redOutliers = outlierData[-which(outlierData$minimalSet == "2088"),]
noOutliers.dge = data.dge[,-which(colnames(data.dge) %in% redOutliers$minimalSet)]

metadataMatching = metadataFiltered[-which(!metadataFiltered$Px_Code %in% colnames(noOutliers.dge)),]
metadataMatching = metadataMatching[match(colnames(noOutliers.dge), metadataMatching$Px_Code),]
head(metadataMatching$Px_Code)
head(colnames(noOutliers.dge))
dim(noOutliers.dge)

palx20.dge = filterByFractionPresent(noOutliers.dge,fraction=0.1,minCount=2)
dim(palx20.dge)
palx20.dge = calcNormFactors(palx20.dge)

minimalSetNoSex = c("Age_mos.","PCT_CORRECT_STRAND_READS","PCT_INTRONIC_BASES", "Initial_growth_duration_days")
tempModel = model.matrix(as.formula(paste("~",paste(minimalSetNoSex,collapse = "+"),sep = "")), data = metadataMatching)
data.voom = voom(palx20.dge,tempModel,plot=FALSE) 
fit = lmFit(data.voom,tempModel)
outData = residuals(fit,y = data.voom)

piezoData  = cbind(outData[which(rownames(outData) == "ENSG00000103335"),], gsub(pattern = "X",replacement = "",x = colnames(outData)))
colnames(piezoData) = c("ENSG00000103335", "Px_Code")
piezoData = as.data.frame(piezoData)
piezoData$ENSG00000103335 = as.numeric(as.character(piezoData$ENSG00000103335))
head(piezoData)

shortMeta = metadataFiltered[-which(is.na(metadataFiltered$SAMPLE)),]

PiezoMeta = merge(shortMeta, allSamples, by.x = "SAMPLE", by.y = "SAMPLE", all.x = FALSE, all.y = FALSE)
head(PiezoMeta)
duplicated(PiezoMeta$SAMPLE)
duplicated(PiezoMeta$Px_Code)

PiezoMerged = merge(PiezoMeta, piezoData, by.x = "Px_Code", by.y = "Px_Code", all.x = FALSE, all.y = FALSE)
head(PiezoMerged)
duplicated(PiezoMerged$SAMPLE)
duplicated(PiezoMerged$Px_Code)
beeswarm(PiezoMerged$ENSG00000103335~PiezoMerged$PIEZOgroup, ylab = "Piezo1 residualized exp")


# group 1 vs control
short = PiezoMerged[-which(PiezoMerged$PIEZOgroup %in% c("2","3")),]
beeswarm(short$ENSG00000103335~short$caseStatus, ylab = "Piezo1 residualized exp", main = "group 1")
t.test(short$ENSG00000103335~short$caseStatus, ylab = "Piezo1 residualized exp")

# groups 1 & 2 vs control
short = PiezoMerged[-which(PiezoMerged$PIEZOgroup == "3"),]
beeswarm(short$ENSG00000103335~short$caseStatus, ylab = "Piezo1 residualized exp", main = "groups 1 & 2")
t.test(short$ENSG00000103335~short$caseStatus, ylab = "Piezo1 residualized exp")


# all Piezo1 variants vs control
boxplot(PiezoMerged$ENSG00000103335~PiezoMerged$caseStatus, main = "all groups")
t.test(PiezoMerged$ENSG00000103335~PiezoMerged$caseStatus)


