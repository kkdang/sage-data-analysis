datafiles = dir()
datafiles_nosummary = datafiles[grep("summary", datafiles, invert = TRUE)]
summaryfiles = datafiles[grep("summary", datafiles)]

readFC=function(inFile){
  print(inFile)
  data = read.delim(inFile, header = TRUE, skip = 1, row.names =1)
  if (ncol(data) > 5) { return(data[,6]) }
  else {return(rep(NA,nrow(data))) }
}

counts = sapply(datafiles_nosummary,readFC)
data = read.delim(datafiles_nosummary[1], header = TRUE, row.names =1, skip = 1)
rownames(counts) = rownames(data)
write.csv(t(counts), file="featurecounts_readcounts.csv", quote = FALSE)



readFCsummary=function(inFile){
  print(inFile)
  data = read.delim(inFile, header = TRUE, row.names =1)
  if (ncol(data) > 0) { return(data[,1])  }
  else {return(rep(NA,9)) }
}

metrics = sapply(summaryfiles,readFCsummary)
data = read.delim(summaryfiles[1], header = TRUE, row.names =1)
rownames(metrics) = rownames(data)
write.csv(t(metrics), file="featurecounts_metrics.csv", quote = FALSE)
# To do a CSV upload, will need to first delete the header line of this file and all rows containing NA




