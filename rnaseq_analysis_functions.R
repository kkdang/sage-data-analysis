# KKD for Sage Bionetworks
# Jan. 14, 2014
# RNAseq analysis functions

library(limma)
library(FactoMineR)


countDetectedByBiotype=function(sampleCounts,biotypes,filter=2){
  tofilter = which(sampleCounts < filter)
  tofilter = c(tofilter, which(sampleCounts == NA))
  filteredCounts = as.matrix(table(biotypes[-tofilter]))
  bioCounts = matrix(0, nrow = length(table(biotypes)), ncol = 1)
  rownames(bioCounts) = rownames(table(biotypes))
  bioCounts[which(rownames(bioCounts) %in% rownames(filteredCounts)),1] = filteredCounts[na.omit(match(rownames(bioCounts), rownames(filteredCounts))),1]
  return(bioCounts[,1])
}

regressDetectionOnMetadata=function(inDetectCounts, inMetaData){
  lmStats = matrix(0, nrow = nrow(inDetectCounts), ncol = 2)
  for (k in 1:nrow(inDetectCounts)) {
    results = summary(lm(inDetectCounts[k,] ~ inMetaData))
    lmStats[k,1] = as.vector(results$coefficients[2,1])
    lmStats[k,2] = as.vector(results$coefficients[2,4])
  }
  return(lmStats)
}

countDetected=function(x, filter=2){ sum(x > filter) }

filterByFractionPresent=function(inCounts,fraction=0.1,minCount=2){
  if (typeof(inCounts) == "DGEList"){
    retain = apply(getCounts(inDGE),MARGIN=1,FUN=palx, x=fraction, filter=minCount)
    filtered = DGEList(getCounts(inDGE[which(retain),]), group = inDGE$samples$group)
  }
  else {
    retain = apply(inCounts,MARGIN=1,FUN=palx, x=fraction, filter=minCount)
    filtered = inCounts[which(retain),]
  }
  return(filtered)
}

palx=function(counts,x=0.1,filter=2){
  tmp = length(which(counts >= filter))/length(counts)
  if (tmp >= x) { return(TRUE)}
  else { return(FALSE)}
}

addDetected=function(counts,filter=2){
  detected = apply(counts,MARGIN=1,function(x) {sum(x > filter)})
  counts$detected = detected
  return(counts)
}

plotExploreData=function(inData,plotTitle){
  if (length(which(!is.na(inData))) == 0) {next}
  if (class(inData) %in% c("character", "logical")) {
    barplot(table(inData), main = plotTitle, las = 2, col = "cornsilk1")
  }
  if (class(inData) %in% c("numeric", "integer")) {
    hist(inData, main=plotTitle, col = "cornsilk1", xlab = "")
  }
  if (class(inData) %in% c("data.frame")) {
    hist(as.numeric(inData), main=plotTitle, col = "cornsilk1", xlab = "")
  }
}
runGoseq=function(lrt,pcutoff=0.05){
  result = topTags(lrt,n=10000)
  plotSmear(lrt,de.tags=rownames(result$table[which(result$table$FDR<pcutoff),]))
  genes=as.integer(p.adjust(result$table$PValue[result$table$logFC!=0], method="BH")<pcutoff)
  names(genes)=row.names(result$table[result$table$logFC!=0,])
  table(genes)
  pwf=nullp(genes,"hg19","ensGene")
  GO.wall=goseq(pwf,"hg19","ensGene")
  enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<pcutoff]
  for(go in enriched.GO){
    x = GOTERM[[go]]
    print(Term(x))
  }
}

calculateDE=function(in.fit,inContrast,plotTitle=modelName,pcutoff=0.05,goseq=FALSE,limma=FALSE,inCoef){
  if (limma){
    fit = eBayes(in.fit)
    if (max(inContrast) > 0) {
      fitContrast = eBayes(contrasts.fit(fit, inContrast))
      result = topTable(fitContrast,number = 10000,sort.by = "p")
      x = length(which(result$adj.P.Val<pcutoff))
      return(list(numSig=x,fit=fitContrast))
    }
    else { 
      result = topTable(fit=fit,coef=inCoef,number = 10000,sort.by = "p")   
      x = length(which(result$adj.P.Val<pcutoff))
      return(list(numSig=x,fit=fit))
    }
  }
  else {
    lrt = glmLRT(in.fit, contrast=inContrast)
    result = topTags(lrt,n=10000)
    x = length(which(result$table$FDR<pcutoff))
    plotSmear(lrt,de.tags=rownames(result$table[which(result$table$FDR<pcutoff),]),main=plotTitle)
    if (goseq) {runGoseq(lrt,pcutoff=pcutoff)}
    return(list(numSig=x,lrt=lrt))
  }
}

plotFCfractionalResults=function(inFC){
  totals = rowSums(inFC)
  fractions = inFC / totals
  op = par(mai = c(3,1.3,1,1))
  boxplot(fractions, las = 2, col = "powderblue", ylab = "fraction of total reads",main = "read destinies")
  par(op)
  op = par(mfrow = c(1,2))
  hist(fractions[,1], main = "distribution of assigned counts", xlab = "fraction of assigned pairs")
  hist(inFC[,1], main = "distribution of assigned counts", xlab = "number of assigned pairs")
  par(op)
}

plotSFfractionalResults=function(inSF){
  fractions = inSF[,c(2,3)] / inSF[,1]
  op = par(mai = c(3,1.3,1,1))
  boxplot(fractions, las = 2, col = "powderblue", ylab = "fraction of total kmers",main = "kmer destinies")
  par(op)
  op = par(mfrow = c(1,2))
  hist(inSF[,4], main = "distribution of mapped kmers", xlab = "fraction of mapped kmers (mappedRatio)")
  hist(inSF[,2], main = "distribution of mapped kmers", xlab = "number of mapped kmers")
  par(op)
}

FM_PCA=function(in.data,plot=TRUE){
  result.prc = PCA(X = t(in.data),scale.unit = TRUE,ncp = 5,graph = FALSE)
  if (plot){ barplot(result.prc$eig$"percentage of variance"[1:20], ylab = "% variance") } 
  return(result.prc)
}

prcomp_PCA=function(in.data,plot=TRUE){
  prc = prcomp(t(in.data),center=TRUE,scale=TRUE,retx = TRUE)
  if (plot) { barplot((prc$sdev^2/sum(prc$sdev^2))[1:20], ylab = "fraction variance") } 
  return(prc)
}

voom_on_mean=function(in.dge,plot=TRUE){
  voom_data = voom(in.dge,plot=plot)
  tmp_fit = lmFit(voom_data,design=rep(1,ncol(voom_data)))
  residuals(tmp_fit,y = voom_data)
}
