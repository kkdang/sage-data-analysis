# Combine Lilly test htseq count data

setwd("~/Computing/NIA-AMP-AD/reprocessed/assay_rosmap_counts/")
x = dir()

samples = unique(sapply(as.list(x),function(x){paste(unlist(strsplit(x,"_"))[1:2],collapse = "_")}))


readHTS=function(inFile){
#  print(inFile)
  data = read.delim(inFile, header = FALSE, row.names =1)
  if (ncol(data) > 0) { return(data[,1]) }
  else {return(rep(NA,nrow(data))) }
}

sumAssayLevel=function(inSample,PEorSE="paired"){
    print(inSample)
    sample_assays = x[grep(pattern = inSample,x = x)]
    samplePairedAssays = sample_assays[grep(PEorSE,sample_assays)]
    countsOriginal = sapply(as.list(samplePairedAssays), readHTS)
    return(rowSums(countsOriginal))
}


pairedSums = sapply(as.list(samples),sumAssayLevel)
singleSums = sapply(as.list(samples),sumAssayLevel,PEorSE="single")

tail(pairedSums)

data = read.delim(x[1], header = FALSE, row.names =1)
rownames(pairedSums) = rownames(data)
rownames(singleSums) = rownames(pairedSums)
colnames(pairedSums) = samples
colnames(singleSums) = samples
head(pairedSums)
head(singleSums)
write.table(pairedSums,file = "../lilly_test_assaylevel_counts.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(singleSums,file = "../lilly_test_assaylevel_countsSE.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


compDataPath = '~/Computing/NIA-AMP-AD/reprocessed/AWS_htseq_counts_test/'
y = dir(compDataPath)
yy = y[grep(pattern = "single",x = y,invert = TRUE)]

compData = sapply(as.list(paste(compDataPath,yy,sep = "")),readHTS)
colnames(compData) = sapply(as.list(yy),function(z){unlist(strsplit(x = z,split = "[.]"))[1]})
tmp = read.delim(paste(compDataPath,yy[1],sep = ""), header = FALSE, row.names =1)
rownames(compData) = rownames(tmp)
head(compData)


zz = y[grep(pattern = "single",x = y,invert = FALSE)]

compDataSE = sapply(as.list(paste(compDataPath,zz,sep = "")),readHTS)
colnames(compDataSE) = sapply(as.list(zz),function(z){unlist(strsplit(x = z,split = "[.]"))[1]})
rownames(compDataSE) = rownames(tmp)
head(compDataSE)




# compare
compDataStats = compData[(nrow(compData)-4):nrow(compData),]
compData = compData[1:(nrow(compData)-5),]
compData = compData[,order(colnames(compData))]
compData = compData[order(rownames(compData)),]


pairedSumsStats = pairedSums[(nrow(pairedSums)-4):nrow(pairedSums),]
pairedSums = pairedSums[1:(nrow(pairedSums)-5),]
pairedSums = pairedSums[,order(colnames(pairedSums))]
pairedSums = pairedSums[order(rownames(pairedSums)),]


head(compData)
head(pairedSums)

which(!colnames(pairedSums) %in% colnames(compData))
pairedSumTrunc = pairedSums[,-1]
tail(compData)
tail(pairedSumTrunc)



diffs = pairedSumTrunc-compData
boxplot(diffs, las = 2)
boxplot(diffs, las = 2, ylim = range(1e4,-1e4))


colsToCompare = match(colnames(compDataSE),colnames(singleSums))
diffsSE = singleSums[,colsToCompare] - compDataSE 
range(diffsSE)
