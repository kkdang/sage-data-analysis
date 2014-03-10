# KKD for Sage Bionetworks
# Oct. 1, 2013
# Retrieves discord matrix file and functions for processing it.

library('synapseClient')
setwd('~/Computing/commonmind/analysis/discord')

synapseLogin()
discordEntity = synGet('syn2342434', version=NULL, downloadFile=T, downloadLocation=getwd(), ifcollision="overwrite.local", load=F)
data = read.delim(getFileLocation(discordEntity))


## functions
convertData=function(inVal){
  return(eval(parse(text=as.character(inVal))))
}

getDenom=function(inVal){
  return(as.numeric(unlist(strsplit(as.character(inVal), "/"))[2]))
}

convertMatrix=function(inData){
  numeric_data = matrix(NA,nrow=nrow(inData),ncol=ncol(inData))
  for (i in 1:NROW(inData)) {
    numeric_data[i,] = unlist(lapply(inData[i,], convertData))
  }
  rownames(numeric_data) = colnames(inData)
  colnames(numeric_data) = colnames(inData)
  hist(numeric_data, xlab = "discordance")
  return(numeric_data)
}


getDenomMatrix=function(inData){
  numeric_data = matrix(NA,nrow=nrow(inData),ncol=ncol(inData))
  for (i in 1:NROW(inData)) {
    numeric_data[i,] = unlist(lapply(inData[i,], getDenom))
  }
  rownames(numeric_data) = colnames(inData)
  colnames(numeric_data) = colnames(inData)
  #  hist(numeric_data)
  return(numeric_data)
}

stripNames=function(inName){
  elements = strsplit(inName,split="_")
  if (length(elements[[1]]) == 4) {
    return(paste(elements[[1]][1], elements[[1]][[4]], sep = "_"))
  }
  else if (length(elements[[1]]) == 5) {
    return(paste(elements[[1]][1], elements[[1]][[3]], elements[[1]][[5]], sep = "_"))
  }
}
## end functions
