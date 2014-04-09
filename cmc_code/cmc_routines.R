### Common CMC routines
  
# Load synapse and Ensembl gene counts file, exclude bad samples

library('synapseClient')
synapseLogin()
dataFile = synGet('syn2340130')
metaFile = synGet('syn2299154')
geneCounts  = read.delim(dataFile@filePath, header=T, row.names = 1)
metadata = read.csv(getFileLocation(metaFile))

#Which samples are repeated in the input file?
which(table(metadata$DLPFC_RNA_isolation..Sample.RNA.ID) > 1)

#Make final dataset, with bad samples excluded. 
metadata_freeze = metadata[-which(metadata$DLPFC_RNA_report..Exclude. == 1),]
total_samples = nrow(metadata_freeze)
table(metadata_freeze$DLPFC_RNA_report..Exclude.)
metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID = factor(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)

metadata_freeze = metadata_freeze[order(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID),]
geneCounts_freeze = geneCounts[,which(as.factor(colnames(geneCounts)) %in% metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)]
geneCounts_freeze = geneCounts_freeze[,order(as.factor(colnames(geneCounts_freeze)))]

head(colnames(geneCounts_freeze))
head(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)
tail(colnames(geneCounts_freeze))
tail(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)
