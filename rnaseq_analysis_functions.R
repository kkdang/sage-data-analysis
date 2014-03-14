# KKD for Sage Bionetworks
# Jan. 14, 2014
# RNAseq analysis functions

library(biomaRt)
#Hs = useMart("ensembl") # use this one normally
Hs=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Hs = useDataset("hsapiens_gene_ensembl", Hs)


getByBiotype=function(biotype="protein_coding", gene=TRUE){
  if (gene==TRUE){
    return(getBM(attributes=c("ensembl_gene_id"),filters="biotype",values=biotype, mart=Hs))    
  }
  else {
    ### need to fix this part
#   return(getBM(attributes=c("ensembl_transcript_id"),filters="biotype",values=biotype, mart=Hs))
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

countDetected=function(x, filter=2){ length(which(x > filter)) }
