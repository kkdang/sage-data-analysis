# KKD for Sage Bionetworks
# Jan. 14, 2014
# RNAseq analysis functions

library(biomaRt)
library(limma)
library(FactoMineR)
#Hs = useMart("ensembl") # use this one normally
Hs=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Hs = useDataset("hsapiens_gene_ensembl", Hs)


getByBiotype=function(biotype="protein_coding", gene=TRUE){
  if (gene==TRUE){
    return(getBM(attributes=c("ensembl_gene_id"),filters="biotype",values=biotype, mart=Hs))    
  }
  else {
   return(getBM(attributes=c("ensembl_transcript_id"),filters="transcript_biotype",values=biotype, mart=Hs))
  }
}

getGeneLengths=function(){
  temp = getBM(attributes=c("ensembl_gene_id", "start_position", "end_position"), mart=Hs) 
  temp$length = temp$end_position - temp$start_position
  return(temp) 
}

getGenesForGOOffspring=function(go_acc,go="CC"){
  library('GO.db')
  temp = switch(go,
                CC = as.list(GOCCOFFSPRING),
                BP = as.list(GOBPOFFSPRING),
                MF = as.list(GOMFOFFSPRING))
  
  offspringGOTerms = temp[[which(names(temp) == go_acc)]]
  allENSG = sapply(offspringGOTerms, getGenesForGOTerm)
  return(unique(unlist(allENSG)))
}


getGenesForGOTerm=function(go_acc){
  return(unique(getBM(attributes=c("ensembl_gene_id"), filters="go_id", values = go_acc, mart=Hs)))
}

densityPlot=function(inDGE,mainLabel="",log=TRUE,ub=1){
  x.cpm = cpm(x=inDGE,normalized.lib.sizes=TRUE,log=log)
  plot(density(x.cpm[,1]), col = "blue", main = mainLabel, ylim = range(0,ub))
  for (i in 2:ncol(x.cpm)) {
    lines(density(x.cpm[,i]), col = "blue")
  }
}

addBiotype=function(inCounts,gene=TRUE){
  if (gene==TRUE){
    biotypes = getBM(attributes=c("ensembl_gene_id", "gene_biotype"),filters="ensembl_gene_id",values=rownames(inCounts), mart=Hs)
  }
  else {
    biotypes = getBM(attributes=c("ensembl_transcript_id", "transcript_biotype"),filters="ensembl_transcript_id",values=rownames(inDGE), mart=Hs)
  }
  inCounts$biotype = rep("NA", nrow(inCounts))
  inCounts$biotype = biotypes$gene_biotype[match(rownames(inCounts), biotypes$ensembl_gene_id)]
  return(inCounts)
}

addEntrez=function(inCounts,gene=TRUE){
  if (gene==TRUE){
    entrez = getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters="ensembl_gene_id",values=rownames(inCounts), mart=Hs)
  }
  else {
    # need to write this condition
#    biotypes = getBM(attributes=c("ensembl_transcript_id", "transcript_biotype"),filters="ensembl_transcript_id",values=rownames(inDGE), mart=Hs)
  }
  inCounts$entrez = rep("NA", nrow(inCounts))
  inCounts$entrez = entrez$entrezgene[match(rownames(inCounts), entrez$ensembl_gene_id)]
  return(inCounts)
}



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

filterByFractionPresent=function(inDGE,fraction=0.1,minCount=2){
  retain = apply(getCounts(inDGE),MARGIN=1,FUN=palx, x=fraction, filter=minCount)
  filtered = DGEList(getCounts(inDGE[which(retain),]), group = inDGE$samples$group)
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

calculateDE=function(in.fit,inContrast,plotTitle=modelName,pcutoff=0.05,goseq=FALSE,limma=FALSE){
  if (limma){
    fit = eBayes(in.fit)
    fitContrast = eBayes(contrasts.fit(fit, inContrast))
    result = topTable(fitContrast,number = 10000,sort.by = "p")
    x = length(which(result$adj.P.Val<pcutoff))
    return(list(numSig=x,fit=fitContrast))
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

FM_PCA=function(in.voom,plot=TRUE){
  result.prc = PCA(X = t(in.voom),scale.unit = TRUE,ncp = 5,graph = FALSE)
  if (plot){ barplot(result.prc$eig$"percentage of variance"[1:20], ylab = "% variance") } 
  return(result.prc)
}

prcomp_PCA=function(in.voom,plot=TRUE){
  prc = prcomp(t(in.voom),center=TRUE,scale=TRUE,retx = TRUE)
  if (plot) { barplot((prc$sdev^2/sum(prc$sdev^2))[1:20], ylab = "fraction variance") } 
  return(prc)
}

voom_normalize=function(in.dge,plot=TRUE){
  voom_data = voom(in.dge,plot=plot)
  tmp_fit = lmFit(voom_data,design=rep(1,ncol(voom_data)))
  residuals(tmp_fit,y = voom_data)
}
