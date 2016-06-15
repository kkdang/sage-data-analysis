# KKD for Sage Bionetworks
# Jan. 14, 2014
# biomart query convenience functions

library(biomaRt)
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

getHGNC=function(inENSG,gene=TRUE){
  geneNames = getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),mart=Hs)    
  if (gene==TRUE){
    return(geneNames[match(inENSG,geneNames[,1]),])    
  }
  else { # write this condition
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


addBiotype=function(inCounts,gene=TRUE){
  if (gene==TRUE){
    biotypes = getBM(attributes=c("ensembl_gene_id", "gene_biotype"),filters="ensembl_gene_id",values=rownames(inCounts), mart=Hs)
  }
  else {
    biotypes = getBM(attributes=c("ensembl_transcript_id", "transcript_biotype"),filters="ensembl_transcript_id",values=rownames(inCounts), mart=Hs)
  }
  inCounts$biotype = rep("NA", nrow(inCounts))
  inCounts$biotype = biotypes$gene_biotype[match(rownames(inCounts), biotypes$ensembl_gene_id)]
  return(inCounts)
}

addEntrez=function(inCounts,gene=TRUE){
  if (gene==TRUE){
    entrez = getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters="ensembl_gene_id",values=rownames(inCounts), mart=Hs)
    inCounts$entrez = rep("NA", nrow(inCounts))
    inCounts$entrez = entrez$entrezgene[match(rownames(inCounts), entrez$ensembl_gene_id)]
  }
  else {
    entrez = getBM(attributes=c("ensembl_transcript_id", "entrezgene"),filters="ensembl_transcript_id",values=rownames(inCounts), mart=Hs)
    inCounts$entrez = rep("NA", nrow(inCounts))
    inCounts$entrez = entrez$entrezgene[match(rownames(inCounts), entrez$ensembl_transcript_id)]
  }
  return(inCounts)
}
