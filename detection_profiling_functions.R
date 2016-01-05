# KKD for Sage Bionetworks
# Jan. 22, 2015
# RNAseq detection profiling functions

library('rGithubClient')
library(edgeR)
sageCode = getRepo(repository="kkdang/sage-data-analysis")
sourceRepoFile(sageCode, "rnaseq_analysis_functions.R")


## Detection analysis of biotypes

### Consolidating and removing some Ensembl biotypes
simplifyBiotypeCategories=function(inData,biomart=Hs){
  # input "inData" is DGEList or raw counts matrix
  if (class(inData) == "DGEList"){
    paloCounts = addBiotype(as.data.frame(getCounts(inData)))
  }
  else { paloCounts = addBiotype(as.data.frame(inData)) }
  origBiotypeCounts = table(getBM(attributes=c("ensembl_gene_id", "gene_biotype"), mart=biomart)[,2])
  # Consolidate total counts of all pseudogene categories into a single category
  allpseudo = grep(pattern = "pseudogene",x = names(origBiotypeCounts))
  tmp = origBiotypeCounts[-allpseudo]
  pseudogene = sum(origBiotypeCounts[allpseudo])
  totalBiotypeCounts = tmp
  totalBiotypeCounts[length(tmp)+1] = pseudogene
  names(totalBiotypeCounts)[length(tmp)+1] = "pseudogene"
  
  # Get rid of observed biotypes with very low counts
  observedBioCounts = table(paloCounts$biotype)
  typesToRemove = c(names(observedBioCounts[which(observedBioCounts < 50)]), "misc_RNA")
  toRemove = which(paloCounts$biotype %in% typesToRemove)
  tmp = paloCounts[-toRemove,]
  
  # Consolidate remaining pseudogene categories into a single category. This undercounts observed pseudogenes by about max(100), since some categories are previously removed for low count. 
  paloCounts_BTFilter = tmp
  allpseudo = grep(pattern = "pseudogene",x = tmp$biotype)
  paloCounts_BTFilter$biotype = replace(tmp$biotype, list = allpseudo, "pseudogene")
  print(table(paloCounts_BTFilter$biotype))
  observedBioCounts = table(paloCounts_BTFilter$biotype)
  return(list(data=paloCounts_BTFilter, totalCounts=totalBiotypeCounts, observedCounts=observedBioCounts))
}

### Plot distribution of number of samples detecting a gene.
plotNumSamplesDetecting=function(inCounts,minCount=1){ # Minimum number of reads to be deemed detected
  paloCounts_wdetect = addDetected(inCounts[,1:ncol(inCounts)-1],filter = minCount-1)
  limits = range(paloCounts_wdetect$detected, na.rm = TRUE)
  hist(paloCounts_wdetect$detected, main = paste("gene detection -- all genes\nmin count =", minCount, "read", sep = " "), xlab=paste("number of samples detecting (", limits[1], "-", limits[2], ")", sep = " "), ylab = "number of genes", col = "lightcyan2")
  paloCounts_wdetect$biotype = inCounts$biotype
  typesToPlot = names(table(paloCounts_wdetect$biotype))
  op = par(mfrow=c(3,3))
  for (i in 1:length(typesToPlot)) {
    hist(paloCounts_wdetect$detected[paloCounts_wdetect$biotype == typesToPlot[i]], main = typesToPlot[i], xlab = "number samples detecting", col = "lightcyan2", xlim = range(0,ncol(inCounts)+1))
  }
  par(op)
  return(paloCounts_wdetect)
}

plotDetectionByBiotype=function(inBiotypeCounts,inNormCounts,inLibSize){
  op = par(mfrow = c(3,3))
  for (i in 1:nrow(inBiotypeCounts)) {
    plot(log10(inLibSize), inNormCounts[i,], main = rownames(inBiotypeCounts)[i], xlab = "log10(mapped reads)", ylab = "fraction detected", ylim = range(0,1))
    lines(lowess(log10(inLibSize), inNormCounts[i,]), col = "blue", lwd = 2)    
  }
  par(op)
}


### Regression analysis of relationship between detection and depth per biotype. Print linear regression results.
runRegressionAndPlot=function(inCounts,bioCategoryCounts,totalBiotypeCounts,inLibSize,minCount=1){
  lmResults = matrix(0, nrow = length(bioCategoryCounts), ncol = 2)
  rownames(lmResults) = names(bioCategoryCounts)
  colnames(lmResults) = c("estimate", "pval")
  
  lmRanges = matrix(0, nrow = length(bioCategoryCounts), ncol = 2)
  rownames(lmRanges) = names(bioCategoryCounts)
  colnames(lmRanges) = c("low", "high")
  
  lmMeds = matrix(0, nrow = length(bioCategoryCounts), ncol = 1)
  rownames(lmRanges) = names(bioCategoryCounts)
  
  j = 1
  total_samples = ncol(inCounts) - 2
  biotypeCountsPerSampleFilter = apply(inCounts[,1:total_samples], 2, countDetectedByBiotype, biotypes = inCounts$biotype, filter=minCount)
  EnsemblTotal = totalBiotypeCounts[match(rownames(biotypeCountsPerSampleFilter), names(totalBiotypeCounts))]
  normBiotypeCounts = apply(biotypeCountsPerSampleFilter, MARGIN=2, function(x) x/EnsemblTotal)
  plotDetectionByBiotype(inBiotypeCounts=biotypeCountsPerSampleFilter[,1:total_samples],inNormCounts=normBiotypeCounts[,1:total_samples],inLibSize = inLibSize)
  coords = match(rownames(normBiotypeCounts), rownames(lmResults))
  lmResults[coords,j:(j+1)] = regressDetectionOnMetadata(inDetectCounts=normBiotypeCounts, inMetaData=log10(inLibSize))
  lmRanges[coords,j:(j+1)] = t(apply(normBiotypeCounts,MARGIN=1,FUN=range))
  lmMeds[coords,1] = apply(normBiotypeCounts,MARGIN=1,FUN=median, na.rm = TRUE)
  j = j + 2
  
  # sets of two columns: first is beta estimate for slope, second is p-value
  print(lmResults)
  print(lmRanges)
  print(lmMeds)
  return(list(estimates=lmResults,range=lmRanges,medians=lmMeds))
}


### Quadrant plot of fraction detected vs regression estimate.
quadrantPlot=function(inResults,observedBioCounts,plotTitle=""){
  par(mfrow = c(1,1))
  plotcols = rainbow(nrow(inResults$estimates))
  plotcex = log10(observedBioCounts)
  indexToPlot = 1 # values = 1,3,5
  altIndex = 1 # values = 1,2,3
  plot(inResults$medians[,altIndex],inResults$estimates[,indexToPlot], xlim = range(0,1), ylim = range(0,0.5), xlab = "fraction detected", ylab = "genes per depth (slope estimate)", col = plotcols, pch = 16,cex=plotcex,main=plotTitle)
  segments(x0=inResults$range[,indexToPlot],y0=inResults$estimates[,indexToPlot],x1=inResults$range[,(indexToPlot+1)],y1=inResults$estimates[,indexToPlot], col = plotcols)
  legend("topright",legend=rownames(inResults$range),col=plotcols,pch=16)
}




## Effect of library quality (eg RIN) on pseudogene detection.

# Setting up a logistic regression, where expression level of parent gene is taken into account. First get information on 'parent' genes of pseudogoenes.
dropSuffix=function(ID){
  unlist(strsplit(as.character(ID),split="[.]"))[1]
}

getPseudoParentData=function(){
  psiFile = synGet('syn3141195')
  psiData  = read.delim(psiFile@filePath,header=TRUE)
  
  psiData_parsedPG = sapply(psiData$Parent.gene,dropSuffix)
  psiData$strippedPG = psiData_parsedPG 
  psiData_parsedPT = sapply(psiData$Parent.transcript,dropSuffix)
  psiData$strippedPT = psiData_parsedPT 
  psiData_parsedID = sapply(psiData$ID,dropSuffix)
  psiData$strippedID = psiData_parsedID 
  
  temp = getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id"),filters="ensembl_transcript_id",values=psiData_parsedID, mart=Hs)
  psiData$geneID = temp$ensembl_gene_id[match(psiData_parsedID, temp$ensembl_transcript_id)]
  return(psiData)
}

# Reduce the input to a set of pseudogenes for which there is an expressed parent gene in the dataset.
getPseudoGenesWithExpParent=function(inData,gene=TRUE){
  if(gene){
    pseudoGenesToConsider = psiData$geneID[which(psiData$strippedPG %in% rownames(inData))]  }
  else{
    pseudoGenesToConsider = psiData$strippedID[which(psiData$strippedPT %in% rownames(inData))]  }
  
  pseudoCounts = inData[which(rownames(inData) %in% pseudoGenesToConsider),]
  
  if(gene){   
    parentsToLocate = psiData$strippedPG[match(rownames(pseudoCounts),psiData$geneID)] }
  else {   
    parentsToLocate = psiData$strippedPT[match(rownames(pseudoCounts),psiData$strippedID)] }
  
  parentCounts = inData[match(parentsToLocate,rownames(inData)) ,]
  return(list(pseudoCounts=pseudoCounts, parentCounts=parentCounts))
}

# results = getPseudoGenesWithExpParent(inData = cpm(data.dge,normalized.lib.sizes = TRUE))
# pseudoCounts = results$pseudoCounts
# parentCounts = results$parentCounts
# rm(results)

# Run the regression and plot results. 
# Send normalized count data into this function
# Requires metadata to be a data fram containing Dx, Mapped.Reads, RIN
runLogistic_wDx_andParent=function(gCounts,parentCounts,metadata){
  F = sapply(gCounts,FUN=detectFilter, simplify=TRUE)
  if ((sum(F) < length(gCounts))&(sum(F) > 0)){
    parentExp = log10(parentCounts)
    parentExp[which(parentExp %in% c("Inf", "-Inf"))] = NA
    fit = try(glm(F~metadata$RIN+log10(metadata$Mapped.Reads)+metadata$Dx+parentExp,family=binomial()),silent=T) 
    if(length(fit) > 1) {
      temp = summary(fit)
      if (nrow(temp$coefficients) < 6) { return(rep(NA,6)) }
      return(c(temp$coefficients[2,4], temp$coefficients[2,1], temp$coefficients[3,4], temp$coefficients[3,1],temp$coefficients[6,4], temp$coefficients[6,1] ))
    }
    else{return(rep(NA, 6))}
  }
  else{return(rep(NA, 6))}
}



#What fraction of pseudogenes tested do not change detection with respect to RIN or parent gene expression? This is an upper bound on valid fraction pseudogenes (ie true observations).
upperBoundPseudo=function(pseudoCounts,parentCounts){
  tpseudocounts = data.frame(t(pseudoCounts))
  tparentcounts = data.frame(t(parentCounts))
  RINeffectByGene_wDx_pseudo = mapply(runLogistic_wDx_andParent, gCounts=tpseudocounts, parentCounts=tparentcounts)
  hist(RINeffectByGene_wDx_pseudo[1,], col = "lightcyan2", main = "pval RIN")
  hist(RINeffectByGene_wDx_pseudo[2,], col = "lightcyan2", main = "est RIN", xlim = range(-100,100), breaks = 2000)
  hist(RINeffectByGene_wDx_pseudo[5,], col = "lightcyan2", main = "pval parent")
  hist(RINeffectByGene_wDx_pseudo[6,], col = "lightcyan2", main = "est parent", xlim = range(-10000,10000), breaks = 200)
  
  length(which(RINeffectByGene_wDx_pseudo[1,] < 0.01)) # RIN
  length(which(RINeffectByGene_wDx_pseudo[3,] < 0.01)) # library size
  length(which(RINeffectByGene_wDx_pseudo[5,] < 0.01)) # parent expression
  
  parentAndRIN = union(which(RINeffectByGene_wDx_pseudo[1,] < 0.01), which(RINeffectByGene_wDx_pseudo[5,] < 0.01)) # either RIN or parent
  length(parentAndRIN)
  nonZeroResult = ncol(RINeffectByGene_wDx_pseudo) - length(which(is.na(RINeffectByGene_wDx_pseudo[1,])==TRUE)) # number of genes for which we were able to test (ie not an NA result)
  nonZeroResult
  fractionPotentiallyValid = 1-(length(parentAndRIN)/nonZeroResult)
  return(fractionPotentiallyValid)
}