x = dir('QUANTsailfish')
y = grep("quant.sf", x)
dataTPM = matrix(nrow=119207,ncol=length(y))
dataRPKM = matrix(nrow=119207,ncol=length(y))
dataKPKM = matrix(nrow=119207,ncol=length(y))
dataEstKmers = matrix(nrow=119207,ncol=length(y))
dataEstReads = matrix(nrow=119207,ncol=length(y))
for (i in 1:length(y)){
  data = read.delim(paste("QUANTsailfish",x[y[i]],sep="/"),skip=4)
  dataTPM[,i] = data[,3]
  dataRPKM[,i] = data[,4]
  dataKPKM[,i] = data[,5]
  dataEstKmers[,i] = data[,6]
  dataEstReads[,i] = data[,7]
}


z = as.list(data[,1])
a = sapply(z,function(x){unlist(strsplit(as.character(x), split = "[|]"))[1]})

zz = as.list(x[y])
b = sapply(zz,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})

colnames(dataTPM) = b
write.csv(x=dataTPM,file="dataTPM.csv",quote=FALSE,row.names = a,col.names=b)
colnames(dataKPKM) = b
write.csv(x=dataKPKM,file="dataKPKM.csv",quote=FALSE,row.names = a)
colnames(dataEstReads) = b
write.csv(x=dataEstReads,file="dataEstReads.csv",quote=FALSE,row.names = a)

