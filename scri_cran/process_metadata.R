# KKD for Sage Bionetworks
# March 2, 2015
# Standard processing of cranio metadata

library(corrgram)

metadataEnt = synGet('syn3241418')
metadata = read.delim(getFileLocation(metadataEnt),header = TRUE)


# Change to monthly batches
data = as.matrix(table(substr(metadata$for_RNA_date_set_up,start = 1,7)))
barplot(t(data),main = "for RNA date set up")
x = as.factor(substr(metadata$for_RNA_date_set_up,start = 1,7))
metadata$for_RNA_date_set_up = x

data = as.matrix(table(substr(metadata$for_RNA_date_plated,start = 1,7)))
barplot(t(data),main = "for RNA date plated")
x = as.factor(substr(metadata$for_RNA_date_plated,start = 1,7))
metadata$for_RNA_date_plated = x

data = as.matrix(table(substr(metadata$for_RNA_date_harvest,start = 1,7)))
barplot(t(data),main = "for RNA date harvest")
x = as.factor(substr(metadata$for_RNA_date_harvest,start = 1,7))
metadata$for_RNA_date_harvest = x


# Change to yearly batches
data = as.matrix(table(substr(metadata$Initial_date_set_up,start = 1,4)))
barplot(t(data),main = "Initial date set up")
x = as.factor(substr(metadata$Initial_date_set_up,start = 1,7))
metadata$Initial_date_set_up = x

data = as.matrix(table(substr(metadata$Initial_date_freeze,start = 1,4)))
barplot(t(data),main = "Initial date freeze")
x = as.factor(substr(metadata$Initial_date_freeze,start = 1,7))
metadata$Initial_date_freeze = x



toFilter = c('File', 'Biomarker_ID', 'Time_Trial', 'LIBRARY', 'PCT_RIBOSOMAL_BASES','READ_GROUP','RIBOSOMAL_BASES','STRAND_SPECIFICITY','SeqSampleName')
metadataFiltered = metadata[,-which(colnames(metadata)%in%toFilter)]
colnames(metadataFiltered)

metadataFiltered[,c(2,7,12)] = log(metadataFiltered[,c(2,7,12)])
metadataFiltered$SAMPLE = as.character(metadataFiltered$SAMPLE)

op = par(mfrow = c(3,3))
for (i in 1:ncol(metadataFiltered)){ 
  if (is.numeric(metadataFiltered[,i])){
    hist(metadataFiltered[,i], main = colnames(metadataFiltered)[i])
  }
}
par(op)
corrgram(metadataFiltered,lower.panel = panel.pts,upper.panel = panel.pie)

rm(data,x,i,op,toFilter)
