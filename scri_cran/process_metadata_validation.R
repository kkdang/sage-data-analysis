# KKD for Sage Bionetworks
# Oct. 31, 2016
# Standard processing of cranio metadata

library(corrgram)

processMetadataVal=function(plot=TRUE){
  metadataTable = synTableQuery('SELECT * FROM syn7486948')
  
  toFilter = c('Sample_ID', 'Order')
  metadataFiltered = metadataTable@values[,-which(colnames(metadataTable@values)%in%toFilter)]
  rownames(metadataFiltered) = metadataTable@values$Sample_ID
  
  metadataFiltered$Plate = as.factor(metadataFiltered$Plate)
  metadataFiltered$Batch = as.factor(metadataFiltered$Batch)
  metadataFiltered$Diagnosis = as.character(metadataFiltered$Diagnosis)
  metadataFiltered$Diagnosis[grep("Coronal",metadataFiltered$Diagnosis)] = "Coronal"
  metadataFiltered$Diagnosis[grep("Lambdoid",metadataFiltered$Diagnosis)] = "Lambdoid"
  metadataFiltered$Diagnosis = as.factor(metadataFiltered$Diagnosis)
  if (plot==TRUE) {
    
    op = par(mfrow = c(3,2))
    hist(metadataFiltered$RNA_Days_to_Harvest)
    hist(log(metadataFiltered$RNA_Days_to_Harvest))
    metadataFiltered$RNA_Days_to_Harvest = log(metadataFiltered$RNA_Days_to_Harvest)
    
    pie(table(metadataFiltered$Plate), main = "Plate")
    pie(table(metadataFiltered$Sex), main = "Sex")
    dotchart(table(metadataFiltered$Batch), main = "Batch")
    
    pie(table(metadataFiltered$Diagnosis), main = "Diagnosis")
    par(op)
    dotchart(table(metadataFiltered$Diagnosis))
    
    pie(table(metadataFiltered$Diagnosis), main = "Diagnosis")
    dotchart(table(metadataFiltered$Diagnosis))
    
    
    corrgram(metadataFiltered,lower.panel = panel.pts,upper.panel = panel.pie)
  }
  return(metadataFiltered)
}
#rm(data,x,i,op,toFilter)
