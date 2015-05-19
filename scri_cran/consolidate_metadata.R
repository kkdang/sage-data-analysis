# KKD for Sage Bionetworks
# Feb. 23, 2015
# Produces consolidated metadata file for this dataset


library('synapseClient')
synapseLogin()
library(scales)
library(lubridate)


fixdate=function(inSplitDate){
  if (as.numeric(inSplitDate[3]) < 15) {inSplitDate[3] = paste("20", inSplitDate[3], sep = "")}
  outDate = paste(inSplitDate[1],inSplitDate[2], inSplitDate[3], sep="/")
  return(outDate)
}
cleanUpDates=function(inDateData){
  splitDates = strsplit(as.character(inDateData),split="/")
  fixedDates = sapply(splitDates,fixdate)
  return(as.Date(mdy(fixedDates)))
}


###############
# Metadata File
##############

# Alignment results
alignEntity = synGet('syn2820392')
align = read.csv(getFileLocation(alignEntity))
colnames(align)

# Clinical data
lookupEnt = synGet("syn2770057")
lookup = read.csv(getFileLocation(lookupEnt))
clinicalTable = synTableQuery('SELECT * FROM syn2823605')
clinicalData = clinicalTable@values
controlSampleIDs = clinicalTable@values$SeqSampleName[which(clinicalTable@values$"Sample Type" == "Control")]

# Clean up hideous column names
colnames(clinicalData) = sapply(strsplit(colnames(clinicalData),split = " "), function(x){paste(x,sep = "",collapse="_")})
colnames(clinicalData) = sapply(strsplit(colnames(clinicalData),split = "[()]"), function(x){paste(x,sep = "",collapse="")})
colnames(clinicalData)

noMatch = which(is.na(match(lookup$UDF.Investigator.Sample.Name, clinicalData$Px_Code)))
lookup[noMatch,]
lookupRev = lookup
lookupRev$UDF.Investigator.Sample.Name = as.character(lookup$UDF.Investigator.Sample.Name)
lookupRev$UDF.Investigator.Sample.Name = replace(x=lookupRev$UDF.Investigator.Sample.Name, list=noMatch, values=c("1071", "1059", "1075"))
noMatch2 = which(is.na(match(lookupRev$UDF.Investigator.Sample.Name, clinicalData$Px_Code)))
rm(noMatch, noMatch2, lookup)

# Make column for case-control status
caseStatus = rep("case", nrow(clinicalData))
caseStatus[which(clinicalData$Sample_Type == "Control")] = "control"
clinicalData$caseStatus = caseStatus

# Clean up date data
temp = cleanUpDates(clinicalData[,8])
clinicalData[,8] = temp
temp = cleanUpDates(clinicalData[,9])
clinicalData[,9] = temp
temp = cleanUpDates(clinicalData[,12])
clinicalData[,12] = temp
temp = cleanUpDates(clinicalData[,13])
clinicalData[,13] = temp
temp = cleanUpDates(clinicalData[,14])
clinicalData[,14] = temp
rm(temp)


tmpMerged = cbind(clinicalData,align[match(clinicalData$SeqSampleName, align$SAMPLE),])
head(tmpMerged)

write.table(tmpMerged,file = "cleanedCombinedMetadata.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
meta = File("cleanedCombinedMetadata.txt",parentId='syn2250990',synapseStore=TRUE)
meta = synStore(meta)

