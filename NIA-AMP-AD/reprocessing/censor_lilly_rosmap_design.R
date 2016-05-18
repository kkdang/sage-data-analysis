# Censor ages in ROSMAP -- Lilly design files

library(synapseClient)
synapseLogin()

setwd("~/Computing/NIA-AMP-AD/reprocessed/")

censorEnt=function(inSynID){
  d1Ent = synGet(inSynID)
  d1 = read.delim(getFileLocation(d1Ent))
  toCensor = grep('age', colnames(d1))
  print(range(d1[,toCensor]))
  valsToCensor = which(d1[,toCensor] >= 90)
  
  d1[valsToCensor,toCensor] = "90+"
  d1[valsToCensor,toCensor]
  print(d1[1:100,toCensor])
  
  outName = paste(strsplit(d1Ent@fileHandle$fileName,split = "[.]")[[1]][1], "censored.txt", sep = "_")
  write.table(d1, file = outName, quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
  synStore(entity = File(path = outName,parentId='syn5810600',synapseStore = TRUE))
}

censorEnt('syn5761336')
censorEnt('syn5761340')
