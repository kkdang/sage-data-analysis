### No bias correction ###
dirName = 'QUANTsailfish'

x = dir(dirName)
y = grep("quant.sf", x)
dataTPM = matrix(nrow=119207,ncol=length(y))
dataEstKmers = matrix(nrow=119207,ncol=length(y))
dataEstReads = matrix(nrow=119207,ncol=length(y))
for (i in 1:length(y)){
  data = read.delim(paste(dirName,x[y[i]],sep="/"),skip=4)
  dataTPM[,i] = data[,3]
  dataEstKmers[,i] = data[,6]
  dataEstReads[,i] = data[,7]
}


z = as.list(data[,1])
a = sapply(z,function(x){unlist(strsplit(as.character(x), split = "[|]"))[1]})

zz = as.list(x[y])
b = sapply(zz,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})

colnames(dataTPM) = b
write.csv(x=dataTPM,file="dataTPM.csv",quote=FALSE,row.names = a,col.names=b)
colnames(dataEstKmers) = b
write.csv(x=dataEstKmers,file="dataEstKmers.csv",quote=FALSE,row.names = a)
colnames(dataEstReads) = b
write.csv(x=dataEstReads,file="dataEstReads.csv",quote=FALSE,row.names = a)



### Bias corrected ###

x = dir(dirName)
y = grep("quant_bias_corrected.sf", x)
dataTPM = matrix(nrow=119207,ncol=length(y))
dataEstKmers = matrix(nrow=119207,ncol=length(y))
dataEstReads = matrix(nrow=119207,ncol=length(y))
for (i in 1:length(y)){
  data = read.delim(paste(dirName,x[y[i]],sep="/"),skip=4)
#  dataTPM[,i] = data[,3]
  dataEstKmers[,i] = data[,6]
#  dataEstReads[,i] = data[,7]
}


z = as.list(data[,1])
a = sapply(z,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})

zz = as.list(x[y])
b = sapply(zz,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})

colnames(dataTPM) = b
write.csv(x=dataTPM,file="dataTPM_biascor.csv",quote=FALSE,row.names = a)
colnames(dataEstKmers) = b
write.csv(x=dataEstKmers,file="dataEstKmers_biascor.csv",quote=FALSE,row.names = a)
colnames(dataEstReads) = b
write.csv(x=dataEstReads,file="dataEstReads_biascor.csv",quote=FALSE,row.names = a)



## metrics files ##
library('synapseClient')
synapseLogin()


allFiles = dir(dirName)
metricsFiles = allFiles[grep("count_info", allFiles)]

readSFsummary=function(inFile){
  print(inFile)
  data = read.delim(paste(dirName,inFile,sep="/"), header = FALSE, row.names =1)
  if (ncol(data) > 0) { return(data[,1])  }
  else {return(rep(NA,4)) }
}

SFmetrics = sapply(metricsFiles,readSFsummary)
data = read.delim(paste(dirName,metricsFiles[1],sep="/"), header = FALSE, row.names =1)
rownames(SFmetrics) = rownames(data)

zz = as.list(colnames(SFmetrics))
b = sapply(zz,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})
colnames(SFmetrics) = b

write.csv(t(SFmetrics), file="sailfish_metrics.csv", quote = FALSE)
# To do a CSV upload, will need to first delete the header line of this file and all rows containing NA
