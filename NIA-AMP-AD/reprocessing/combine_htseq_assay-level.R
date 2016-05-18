#! /usr/bin/env Rscript
# KKD for Sage Bionetworks
# Combine Lilly test htseq count data
# May 17, 2016

library(argparse)

parser = ArgumentParser(description='Collect and sum assay-level counts into a sample-by-gene matrix file.')
parser$add_argument('--outFileName', type="character", required=TRUE, help='Prefix for output files.')
parser$add_argument('--wd', default = getwd(), type="character", help='Directory for output files [default %(default)s].')


args = parser$parse_args()


setwd(args$wd)
x = dir()

samples = unique(sapply(as.list(x),function(x){paste(unlist(strsplit(x,"_"))[1:2],collapse = "_")}))
head(samples)


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
write.table(pairedSums,file = paste(args$outFilePrefix, "counts_matrix.txt", sep = "_"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(singleSums,file = paste(args$outFilePrefix, "counts_matrixSE.txt", sep = "_"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


