# KKD for Sage Bionetworks
# Oct. 31, 2016
# Standard processing of cranio metadata

library(corrgram)

processMetadataVal=function(plot=TRUE){
  metadataTable = as.data.frame(synTableQuery('SELECT * FROM syn8548351',includeRowIdAndRowVersion = FALSE))
  lookup = as.data.frame(synTableQuery('SELECT * FROM syn7477113',includeRowIdAndRowVersion=FALSE))
  picardMetrics = read.delim(synGet('syn8489352')$path,row.names = 1, header = TRUE,colClasses = c("character",rep("numeric",21)))
  
  # join all
  picardMetrics$sampName = sapply(rownames(picardMetrics),function(x) {unlist(strsplit(x,split = "[.]"))[1]})
  combined = merge(lookup,y = picardMetrics,by.x = "Sample_Name",by.y = "sampName")
  
  toFilter = c('Order', 'Space')
  metadataFiltered = metadataTable[,-which(colnames(metadataTable)%in%toFilter)]
  allTogetherNow = merge(x=combined,y=metadataFiltered,by.x = "Investigator_Sample_Name",by.y = 'Sample_ID')
  rm(combined)
  
  b = which(colnames(allTogetherNow) == "MEDIAN_5PRIME_TO_3PRIME_BIAS")
  for (i in 3:b) {
    allTogetherNow[,i] = as.numeric(as.character(allTogetherNow[,i]))}
  
  allTogetherNow$Plate = as.factor(allTogetherNow$Plate)
  allTogetherNow$Batch = as.factor(allTogetherNow$Batch)
  allTogetherNow$Diagnosis = as.character(allTogetherNow$Diagnosis)
  allTogetherNow$Diagnosis[grep("Coronal",allTogetherNow$Diagnosis)] = "Coronal"
  allTogetherNow$Diagnosis[grep("Lambdoid",allTogetherNow$Diagnosis)] = "Lambdoid"
  allTogetherNow$Diagnosis = as.factor(allTogetherNow$Diagnosis)
  allTogetherNow$RNA_days_to_Harvest = as.numeric(as.character(allTogetherNow$RNA_days_to_Harvest))
  if (plot==TRUE) {
    
    op = par(mfrow = c(3,2))
    hist(allTogetherNow$RNA_days_to_Harvest)
    hist(log(allTogetherNow$RNA_days_to_Harvest))
    allTogetherNow$RNA_days_to_Harvest = log(allTogetherNow$RNA_days_to_Harvest)
    
    pie(table(allTogetherNow$Plate), main = "Plate")
    pie(table(allTogetherNow$Sex), main = "Sex")
    dotchart(table(allTogetherNow$Batch), main = "Batch")
    
    pie(table(allTogetherNow$Diagnosis), main = "Diagnosis")
    par(op)
    dotchart(table(allTogetherNow$Diagnosis))
    
    pie(table(allTogetherNow$Diagnosis), main = "Diagnosis")
    dotchart(table(allTogetherNow$Diagnosis))
    
  
    allTogetherNowShort = allTogetherNow[,-which(colnames(allTogetherNow) %in% c("RIBOSOMAL_BASES", "PCT_RIBOSOMAL_BASES"))]
    c = which(colnames(allTogetherNowShort) == "INCORRECT_STRAND_READS")
    for (i in 3:c) {
      allTogetherNowShort[,i] = log(allTogetherNowShort[,i]) }
    corrgram(allTogetherNowShort[,3:ncol(allTogetherNowShort)],lower.panel = panel.conf,upper.panel = panel.shade,type = "data",cor.method = "spearman", order = TRUE,label.srt = 40,cex.labels = 0.4)
  }
  return(allTogetherNow)
}
