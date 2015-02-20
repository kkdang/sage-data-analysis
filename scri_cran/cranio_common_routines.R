# KKD for Sage Bionetworks
# Jan. 22, 2015
# Common cranio data analysis functions


library('synapseClient')
synapseLogin()
library(scales)
library(lubridate)
library(edgeR)
library(RColorBrewer)
library('rGithubClient')
sageCode = getRepo(repository="kkdang/sage-data-analysis")
sourceRepoFile(sageCode, "rnaseq_analysis_functions.R")

# Alignment results
alignEntity = synGet('syn2820392')
align = read.csv(getFileLocation(alignEntity))
colnames(align)


#Sample-level expression values
# Fix the sample names of the estimated read count data
lookupEnt = synGet("syn2770057")
lookup = read.csv(getFileLocation(lookupEnt))
clinicalTable = synTableQuery('SELECT * FROM syn2823605')
clinicalData = clinicalTable@values

colnames(clinicalData) = sapply(strsplit(colnames(clinicalData),split = " "), function(x){paste(x,sep = "",collapse="_")})
colnames(clinicalData) = sapply(strsplit(colnames(clinicalData),split = "[()]"), function(x){paste(x,sep = "",collapse="")})


noMatch = which(is.na(match(lookup$UDF.Investigator.Sample.Name, clinicalData$Px_Code)))
lookup[noMatch,]

lookupRev = lookup
lookupRev$UDF.Investigator.Sample.Name = as.character(lookup$UDF.Investigator.Sample.Name)
lookupRev$UDF.Investigator.Sample.Name = replace(x=lookupRev$UDF.Investigator.Sample.Name, list=noMatch, values=c("1071", "1059", "1075"))
noMatch2 = which(is.na(match(lookupRev$UDF.Investigator.Sample.Name, clinicalData$"Px Code")))
rm(noMatch, noMatch2, lookup)


#Clinical and experimental data
controlSampleIDs = clinicalTable@values$SeqSampleName[which(clinicalTable@values$"Sample Type" == "Control")]

fixdate=function(inSplitDate){
  if (as.numeric(inSplitDate[3]) < 15) {inSplitDate[3] = paste("20", inSplitDate[3], sep = "")}
  outDate = paste(inSplitDate[1],inSplitDate[2], inSplitDate[3], sep="/")
  return(outDate)
}
cleanUpDates=function(inDateData){
  splitDates = strsplit(as.character(inDateData),split="/")
  fixedDates = sapply(splitDates,fixdate)
  return(as.Date(mdy(fixedDates)))
}


#Function to generate data object
generateDataAndClinicalObj=function(synId){
  dataEntity = synGet(synId) # estimated Reads
  library('R.utils')
  gunzip(getFileLocation(dataEntity),overwrite=TRUE)
  x = unlist(strsplit(getFileLocation(dataEntity), split = '.gz'))
  detach("package:R.utils", unload=TRUE)
  estReads = read.csv(x, row.names = 1)
  colnames(estReads) = lookupRev$UDF.Investigator.Sample.Name[match(colnames(estReads), paste("X", lookupRev$Sample.Name,sep = ""))]
  
  data.dge = DGEList(counts=estReads,remove.zeros=TRUE)
  
  
  # Remove clinical data for samples that don't have matching seq data.
  noMatch = which(is.na(match(clinicalData$"Px_Code", colnames(data.dge))))
  clinicalDataR = clinicalData[-noMatch,]
  
  # clean up date data
  x = cleanUpDates(clinicalDataR[,8])
  clinicalDataR[,8] = x
  x = cleanUpDates(clinicalDataR[,9])
  clinicalDataR[,9] = x
  x = cleanUpDates(clinicalDataR[,12])
  clinicalDataR[,12] = x
  x = cleanUpDates(clinicalDataR[,13])
  clinicalDataR[,13] = x
  x = cleanUpDates(clinicalDataR[,14])
  clinicalDataR[,14] = x
  
  return(list(DGE=data.dge, VARS=clinicalDataR))
}


pcaCorrelation=function(x,y){cor.test(x,y)$p.value }
pcaCorrelationVal=function(x,y){cor.test(x,y)$estimate }
pcaLM=function(x,y){summary(lm(x ~ y))$coefficients[2,4]}
pcaLMEst=function(x,y){summary(lm(x ~ y))$coefficients[2,1]}

calc_plot_PC_corr=function(in.pca,inClinical=clinicalDataR,categorical=c(),k=15){
  # Using correlation against eigenvectors from all PALO data
  clinicalPCAcorrelations = matrix(NA, nrow = ncol(in.pca), ncol = ncol(inClinical))
  colnames(clinicalPCAcorrelations) = rep("default", ncol(clinicalPCAcorrelations))
  clinicalPCAcorrelationEst = matrix(NA, nrow = ncol(in.pca), ncol = ncol(inClinical))
  colnames(clinicalPCAcorrelationEst) = rep("default", ncol(clinicalPCAcorrelations))
  
  # Correlate with numeric factors using pearson
  toRun = setdiff(seq(1,ncol(inClinical)),categorical)
  for (j in 1:length(toRun)) {
    clinicalPCAcorrelations[,j] = apply(in.pca,MARGIN=2,pcaCorrelation,y=inClinical[,toRun[j]])
    colnames(clinicalPCAcorrelations)[j] = colnames(inClinical)[toRun[j]]
    clinicalPCAcorrelationEst[,j] = apply(in.pca,MARGIN=2,pcaCorrelationVal,y=inClinical[,toRun[j]])
    colnames(clinicalPCAcorrelationEst)[j] = colnames(inClinical)[toRun[j]]
  }
  
  # Regress against categorical factors
  temp = matrix(NA, nrow = ncol(in.pca), ncol = length(categorical))
  colnames(temp) = rep("default", ncol(temp))
  tempEst = matrix(NA, nrow = ncol(in.pca), ncol = length(categorical))
  colnames(tempEst) = rep("default", ncol(tempEst))
  for (j in 1:length(categorical)) {
    temp[,j] = apply(in.pca,MARGIN=2,pcaLM,y=inClinical[,categorical[j]])
    colnames(temp)[j] = colnames(inClinical)[categorical[j]]
    tempEst[,j] = apply(in.pca,MARGIN=2,pcaLMEst,y=inClinical[,categorical[j]])
    colnames(tempEst)[j] = colnames(inClinical)[categorical[j]]
  }
  clinicalPCAcorrelations[,(ncol(clinicalPCAcorrelations)-length(categorical)+1):ncol(clinicalPCAcorrelations)] = temp
  colnames(clinicalPCAcorrelations)[(ncol(clinicalPCAcorrelations)-length(categorical)+1):ncol(clinicalPCAcorrelations)] = colnames(temp)
  clinicalPCAcorrelationEst[,(ncol(clinicalPCAcorrelationEst)-length(categorical)+1):ncol(clinicalPCAcorrelationEst)] = tempEst
  colnames(clinicalPCAcorrelationEst)[(ncol(clinicalPCAcorrelationEst)-length(categorical)+1):ncol(clinicalPCAcorrelationEst)] = colnames(tempEst)
  
  op = par(mai = c(3,1,1,1))
  #  boxplot(clinicalPCAcorrelations[,1:26], las = 2)
  boxplot(clinicalPCAcorrelations, las = 2)
  par(op)
  #  boxplot(t(clinicalPCAcorrelations), las = 2, main = "Correlation of each PCs vs each factor", ylab = "raw pvalue")
  dotchart(t(clinicalPCAcorrelations[1:2,]), xlab = "pvalue", main = "Correlation with PC1 and PC2",bg="tomato",lcolor="tomato",cex=0.9)
  dotchart(t(clinicalPCAcorrelations[1:2,]), xlab = "pvalue", main = "Correlation with PC1 and PC2", xlim = range(0,0.1),bg="tomato",lcolor="tomato",cex=0.9)
  abline(v=0.05, lty = 2)
  abline(v = 0.01, lty = 3)
  #apply(clinicalPCAcorrelations[1:20,],MARGIN=1,stripchart, las = 2)
  
  hist(clinicalPCAcorrelations)
  reds=brewer.pal(7,"Reds")
  dataFiltered = clinicalPCAcorrelations
  dataFiltered[which(clinicalPCAcorrelations > cutoff)] = NA
  heatmap(t(1-dataFiltered[1:k,]),Rowv = NA,Colv = NA,scale = "none",col=reds,margins = c(3,13))
  rdbu=brewer.pal(8,"RdBu")
  dataFiltered = clinicalPCAcorrelationEst
  dataFiltered[which(clinicalPCAcorrelations > cutoff)] = NA
  heatmap(t(dataFiltered[1:k,1:(ncol(dataFiltered)-5)]),Rowv = NA,Colv = NA,scale = "none",col=rdbu,margins = c(3,13))
}