# Removing samples from full dataset

library('synapseClient')
synapseLogin()
dataFile = synGet('syn2340130')

geneCounts  = read.delim(dataFile@filePath, header=T, row.names = 1)

toExclude = c("MSSM_RNA_PFC_72", "MSSM_RNA_PFC_73", "MSSM_RNA_PFC_74", "MSSM_RNA_PFC_214", "PENN_RNA_PFC_82", "MSSM_RNA_PFC_21", "MSSM_RNA_PFC_27", "MSSM_RNA_PFC_385", "MSSM_RNA_PFC_93")

ncol(geneCounts)
geneCounts_release = geneCounts[,-which(as.factor(colnames(geneCounts)) %in% toExclude)]
ncol(geneCounts_release)

head(colnames(geneCounts_release))
tail(colnames(geneCounts_release))

grep(toExclude,colnames(geneCounts))
grep(toExclude,colnames(geneCounts_release))

outFileName = "geneCounts_release"
write.table(x = geneCounts_release,file = outFileName, quote = FALSE,row.names = TRUE, col.names = TRUE)
library(R.utils)
gzip(outFileName)

fileEntity = File(path=paste(outFileName, "gz", sep = "."), parentId='syn3158410')
synSetAnnotations(fileEntity) = list(normalizationStatus=FALSE, summaryLevel="gene", geneIdentifiers="Ensembl")
fileEntity = synStore(fileEntity)
