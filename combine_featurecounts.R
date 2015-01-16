datafiles = dir('COUNTfeaturecounts')
datafiles_nosummary = datafiles[grep("summary", datafiles, invert = TRUE)]
summaryfiles = datafiles[grep("summary", datafiles)]

readFC=function(inFile){
  print(inFile)
  data = read.delim(paste("COUNTfeaturecounts",inFile,sep="/"), header = TRUE, skip = 1, row.names =1)
  if (ncol(data) > 5) { return(data[,6]) }
  else {return(rep(NA,nrow(data))) }
}

counts = sapply(datafiles_nosummary,readFC)
data = read.delim(paste("COUNTfeaturecounts",datafiles_nosummary[1],sep="/"), header = TRUE, row.names =1, skip=1)
rownames(counts) = rownames(data)

zz = as.list(colnames(counts))
b = sapply(zz,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})
colnames(counts) = b
write.csv(counts, file="featurecounts_readcounts.csv", quote = FALSE)



readFCsummary=function(inFile){
  print(inFile)
  data = read.delim(paste("COUNTfeaturecounts",inFile,sep="/"), header = TRUE, row.names =1)
  if (ncol(data) > 0) { return(data[,1])  }
  else {return(rep(NA,9)) }
}

metrics = sapply(summaryfiles,readFCsummary)
data = read.delim(paste("COUNTfeaturecounts",summaryfiles[1],sep="/"), header = TRUE, row.names =1)
rownames(metrics) = rownames(data)

zz = as.list(colnames(metrics))
b = sapply(zz,function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})
colnames(metrics) = b

write.csv(t(metrics), file="featurecounts_metrics.csv", quote = FALSE)
# To do a CSV upload, will need to first delete the header line of this file and all rows containing NA




