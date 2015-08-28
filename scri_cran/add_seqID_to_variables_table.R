# KKD for Sage Bionetworks
# January 9, 2015
# Add seq ID to clinical/exp data table


library('synapseClient')
synapseLogin()
setwd('/Users/kristen/Computing/cranio/')


lookupEnt = synGet("syn2770057")
lookup = read.csv(getFileLocation(lookupEnt))

clinicalTable = synTableQuery('SELECT * FROM syn2823605')
clinicalData = clinicalTable@values

noMatch = which(is.na(match(lookup$UDF.Investigator.Sample.Name, clinicalData$"Px Code")))
lookup[noMatch,]
lookupRev = lookup
lookupRev$UDF.Investigator.Sample.Name = as.character(lookup$UDF.Investigator.Sample.Name)
lookupRev$UDF.Investigator.Sample.Name = replace(x=lookupRev$UDF.Investigator.Sample.Name, list=noMatch, values=c("1071", "1059", "1075"))
noMatch2 = which(is.na(match(lookupRev$UDF.Investigator.Sample.Name, clinicalData$"Px Code")))
rm(noMatch, noMatch2, lookup, lookupEnt)


head(clinicalTable@values)
clinicalTable@values$SeqSampleName = lookupRev$Sample.Name[match(clinicalTableData@values$"Px Code", lookupRev$UDF.Investigator.Sample.Name)]
class(clinicalTable@values$SeqSampleName) = "character"
clinicalData = synStore(clinicalTable,retrieveData = TRUE)
