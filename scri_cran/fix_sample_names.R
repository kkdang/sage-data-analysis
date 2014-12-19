# ---
# title: "fix_sample_names.Rmd"
# author: "Kristen Dang"
# date: "12/18/2014"
# output: html_document
# ---
  
library('synapseClient')
synapseLogin()


# Fix the sample names of the estimated read count data
dataEntity = synGet('syn2961661',version = 1)
library('R.utils')
gunzip(getFileLocation(dataEntity),overwrite=TRUE)
x = unlist(strsplit(getFileLocation(dataEntity), split = '.gz'))
detach("package:R.utils", unload=TRUE)
estReads = read.csv(x, row.names = 1)

lookupEnt = synGet("syn2770057")
lookup = read.csv(getFileLocation(lookupEnt))
clinicalEntity = synGet('syn2731158')
clinicalData = read.csv(getFileLocation(clinicalEntity))

noMatch = which(is.na(match(lookup$UDF.Investigator.Sample.Name, clinicalData$Px.Code)))
lookup[noMatch,]

lookupRev = lookup
lookupRev$UDF.Investigator.Sample.Name = as.character(lookup$UDF.Investigator.Sample.Name)
lookupRev$UDF.Investigator.Sample.Name = replace(x=lookupRev$UDF.Investigator.Sample.Name, list=noMatch, values=c("1071", "1059", "1075"))
noMatch2 = which(is.na(match(lookupRev$UDF.Investigator.Sample.Name, clinicalData$Px.Code)))
rm(noMatch, noMatch2, lookup, lookupEnt)

colnames(estReads) = lookupRev$UDF.Investigator.Sample.Name[match(colnames(estReads), paste("X", lookupRev$Sample.Name,sep = ""))]


save(estReads, file = getFileLocation(dataEntity), compress="gzip")
dataEntity = synStore(dataEntity)
