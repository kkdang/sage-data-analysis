# KKD for Sage Bionetworks
# Oct. 1, 2013
# Converts discord values to numeric, plots and prints potential duplicates to a file

library('synapseClient')
setwd('~/Computing/commonmind/analysis/discord')

synapseLogin()
discordEntity = synGet('syn2342434', version=NULL, downloadFile=T, downloadLocation=getwd(), ifcollision="overwrite.local", load=F)
data = read.delim(getFileLocation(discordEntity))


## functions
convertData=function(inVal){
  return(eval(parse(text=as.character(inVal))))
}

getDenom=function(inVal){
  return(as.numeric(unlist(strsplit(as.character(inVal), "/"))[2]))
}

convertMatrix=function(inData){
  numeric_data = matrix(NA,nrow=nrow(inData),ncol=ncol(inData))
  for (i in 1:NROW(inData)) {
    numeric_data[i,] = unlist(lapply(inData[i,], convertData))
  }
  rownames(numeric_data) = colnames(inData)
  colnames(numeric_data) = colnames(inData)
  hist(numeric_data)
  return(numeric_data)
}


getDenomMatrix=function(inData){
  numeric_data = matrix(NA,nrow=nrow(inData),ncol=ncol(inData))
  for (i in 1:NROW(inData)) {
    numeric_data[i,] = unlist(lapply(inData[i,], getDenom))
  }
  rownames(numeric_data) = colnames(inData)
  colnames(numeric_data) = colnames(inData)
#  hist(numeric_data)
  return(numeric_data)
}

stripNames=function(inName){
  elements = strsplit(inName,split="_")
  if (length(elements[[1]]) == 4) {
    return(paste(elements[[1]][1], elements[[1]][[4]], sep = "_"))
  }
  else if (length(elements[[1]]) == 5) {
    return(paste(elements[[1]][1], elements[[1]][[3]], elements[[1]][[5]], sep = "_"))
  }
}
## end functions


completeData = convertMatrix(data)

matchingFromFullFile = shortNewData[which(rownames(shortNewData) %in% rownames(discordData)),which(colnames(shortNewData) %in% colnames(discordData))]
matchingFromLimitedFile = discordData[which(rownames(discordData) %in% rownames(matchingFromFullFile)),which(colnames(discordData) %in% colnames(matchingFromFullFile))]


matchingFromFullFile = matchingFromFullFile[,order(colnames(matchingFromFullFile))]
matchingFromLimitedFile = matchingFromLimitedFile[,order(colnames(matchingFromLimitedFile))]

hist(matchingFromFullFile-matchingFromLimitedFile)

# find first RNAseq sample
x = min(grep(pattern="_RNA_",x=colnames(completeData)))
end = ncol(completeData)

# look at rna-rna matches only
rna_data_only = completeData[x:end,x:end]
rna_upper = rna_data_only
rna_upper[lower.tri(rna_upper)] = NA
pdf(file="hist_rna-rna_discord.pdf")
hist(rna_upper, main = "RNA-RNA comparisons", xlab = "discordance rate", col = "moccasin", breaks = 40)
dev.off()
pdf(file="hist_rna-rna_discord_zoom.pdf")
hist(rna_upper, main = "RNA-RNA comparisons", xlab = "discordance rate", col = "moccasin", breaks = 40, xlim = range(0.03,0.3), ylim = range(0,100))
dev.off()
range(rna_upper, na.rm = TRUE)
# get a list containing the matches for each sample
rna_dups = apply(rna_upper, 1, function(x) which((x < 0.15) == TRUE))
# how many matches does each sample have?
num_matches = sapply(rna_dups, length)
# which ones have more than one match?
dup_indicies = which(num_matches > 1)
write('RNAseq-RNAseq sample matches\n', file='duplicate_results_summary.txt',append=F)
for (i in 1:length(dup_indicies)) {
  write(paste(names(rna_dups[dup_indicies[i]]), sapply(rna_dups[dup_indicies[i]],names), sep = "-match-"), file = 'duplicate_results_summary.txt', append = T)
}
for (i in 1:length(dup_indicies)) {
  stripchart(rna_upper[dup_indicies[i],])
}


# look at dna-rna matches 
mixed_upper = completeData
mixed_upper[lower.tri(mixed_upper)] = NA
mixed_data_only = mixed_upper[1:(x-1),x:end]
pdf(file="hist_dna-rna_discord.pdf")
hist(mixed_data_only, main = "DNA-RNA comparisons", xlab = "discordance rate", col = "moccasin", breaks = 40)
dev.off()
pdf(file="hist_dna-rna_discord_zoom.pdf")
hist(mixed_data_only, main = "DNA-RNA comparisons", xlab = "discordance rate", col = "moccasin", breaks = 40, xlim = range(0,0.5), ylim = range(0,300))
dev.off()
mixed_dups = apply(mixed_data_only, 1, function(x) which((x < 0.25) == TRUE))
num_matches = sapply(mixed_dups, length)
dup_indicies = which(num_matches >= 1)
self_matches = c()
nonself_matches = c()

write('\n\nRNAseq-DNA chip sample matches\n', file='duplicate_results_summary.txt',append=T)
for (i in 1:length(dup_indicies)) {
  matches = sapply(names(mixed_dups[[dup_indicies[i]]]),stripNames) # strips "RNA" and "PFC" out of match names so they can be compared to genotype names
  identical = which(matches == toupper(names(mixed_dups[dup_indicies[i]])))
  if (length(identical) > 0) { 
    self_matches = c(self_matches, mixed_data_only[dup_indicies[i], mixed_dups[[dup_indicies[i]]][identical]])
    matches = matches[-identical] # remove matching ones
  }
  if (length(matches) > 0) {
    nonself_matches = c(nonself_matches, mixed_data_only[dup_indicies[i], mixed_dups[[dup_indicies[i]]]])
    write(paste(names(mixed_dups[dup_indicies[i]]), matches, sep = "-match-"), file = 'duplicate_results_summary.txt', append = T)
#  write(paste(names(mixed_dups[dup_indicies[i]]), sapply(mixed_dups[dup_indicies[i]],names), sep = "-match-"), file = 'duplicate_results_summary.txt', append = T)
  }
}

pdf(file="density_discordance_matches.pdf")
plot(density(nonself_matches))
lines(density(self_matches), col = "orange")
legend(x="topright", legend=c("nonself", "self"),col=c("black", "orange"), pch = 16)
dev.off()

pdf(file="hist2x_discordance_matches.pdf")
par(mfrow = c(1,2))
hist(nonself_matches, xlim = range(0,0.15), ylim = range(0,160), col = "gray", breaks = 10)
hist(self_matches, xlim = range(0,0.15), ylim = range(0,160), col = "orange", ylab = "Frequency")
par(mfrow = c(1,1))
dev.off()


#for (i in 1:length(dup_indicies)) {
#  stripchart(mixed_data_only[dup_indicies[i],])
#}
matCols = c(161, 203, 100, 125, 101, 227, 426)
# new matches
indicies = c(23, 176, 152, 178, 144, 322, 88)
samples = c("DNA MSSM_27", "DNA MSSM_260", "DNA MSSM_184", "DNA MSSM_225", "DNA MSSM_183", "DNA MSSM_304", "DNA_MSSM_82")
plotStripChart=function(inIndex,inTitle){
  pdf(file=paste("stripchart_", inTitle, ".pdf", sep = ""))
  stripchart(mixed_data_only[inIndex,], main = inTitle)
  dev.off()
}
for (i in 1:length(indicies)){
  plotStripChart(inIndex=indicies[i],inTitle=samples[i])
}



# look at dna-rna matches - number of matched sites
completeDenom = getDenomMatrix(local)

denom_upper = completeDenom
denom_upper[lower.tri(denom_upper)] = NA
denom_data_only = denom_upper[1:(x-1),x:end]
pdf(file="hist_dna-rna_denom.pdf")
hist(denom_data_only, main = "DNA-RNA comparisons", xlab = "num sites", col = "moccasin", breaks = 40)
dev.off()


self_denom = c()
nonself_denom = c()
for (i in 1:length(dup_indicies)) {
  matches = sapply(names(mixed_dups[[dup_indicies[i]]]),stripNames) # strips "RNA" and "PFC" out of match names so they can be compared to genotype names
  identical = which(matches == toupper(names(mixed_dups[dup_indicies[i]])))
  if (length(identical) > 0) { 
    self_denom = c(self_denom, denom_data_only[dup_indicies[i], mixed_dups[[dup_indicies[i]]][identical]])
    matches = matches[-identical] # remove matching ones
  }
  if (length(matches) > 0) {
    nonself_denom = c(nonself_denom, denom_data_only[dup_indicies[i], mixed_dups[[dup_indicies[i]]]])
  }
}
pdf(file="density_denom_matches.pdf")
plot(density(nonself_denom))
lines(density(self_denom), col = "orange")
legend(x="topright", legend=c("nonself", "self"),col=c("black", "orange"), pch = 16)
dev.off()


pdf(file="hist2x_denom_matches.pdf")
par(mfrow = c(1,2))
hist(nonself_denom, xlim = range(0,1500), ylim = range(0,160), col = "gray", breaks = 7)
hist(self_denom, xlim = range(0,1500), ylim = range(0,160), col = "orange")
par(mfrow = c(1,1))
dev.off()



# new matches
plotStripChart=function(inIndex,inTitle){
  pdf(file=paste("stripchart_denom_", inTitle, ".pdf", sep = ""))
  stripchart(denom_data_only[inIndex,], main = inTitle)
  dev.off()
}
for (i in 1:length(indicies)){
  plotStripChart(inIndex=indicies[i],inTitle=samples[i])
}

# random plots
stripchart(denom_data_only[179,])
stripchart(denom_data_only[180,])
stripchart(denom_data_only[181,])

pdf(file="dotchart_denom_newmatches.pdf")
denomsForNewMatches = c()
for (i in 1:length(indicies)) { denomsForNewMatches = c(denomsForNewMatches, denom_data_only[indicies[i], matCols[i]])}
names(denomsForNewMatches) = rownames(denom_data_only)[indicies]
dotchart(denomsForNewMatches)
dev.off()

# Are there ~ 6 RNA samples that have low calls across all genotypes? Are these the six new ones?
rnaDenomMeans = colMeans(denom_data_only)
hist(rnaDenomMeans)
length(which(rnaDenomMeans < 500))
which(rnaDenomMeans < 500)
# There are 8 samples with total calls < 500. But these are low across all samples, not just the six new ones.

## What is the cause of the left-hand tail on the denominator size? Quality problem with the RNA samples?




#######
# check for self-matches in dna-rna comparison 
mixed_upper = completeData
mixed_upper[lower.tri(mixed_upper)] = NA
mixed_data_only = mixed_upper[1:(x-1),x:end]
mixed_dups = apply(mixed_data_only, 1, function(x) which((x < 0.25) == TRUE))
num_matches = sapply(mixed_dups, length)
dup_indicies = which(num_matches >= 1)
missing_self = c()

#write('\n\nRNAseq-DNA chip samples WITHOUT self-matches\n', file='duplicate_results_summary.txt',append=T)
write('\n\nRNAseq-DNA chip samples WITHOUT self-matches\n', file='testfiles.txt',append=T)
for (i in 1:length(dup_indicies)) {
  matches = sapply(names(mixed_dups[[dup_indicies[i]]]),stripNames) # strips "RNA" and "PFC" out of match names so they can be compared to genotype names
  identical = which(matches == toupper(names(mixed_dups[dup_indicies[i]])))
  if (length(identical) == 0) { 
    missing_self = c(missing_self, names(mixed_dups[dup_indicies[i]]))
    write(paste(names(mixed_dups[dup_indicies[i]]), "NO MATCH", sep = "-match-"), file = 'testfiles.txt', append = T) }
#    write(paste(names(mixed_dups[dup_indicies[i]]), "NO MATCH", sep = "-match-"), file = 'duplicate_results_summary.txt', append = T)  }
}

nodup_indicies = which(num_matches == 0)


for (i in 1:length(nodup_indicies)) {
  missing_self = c(missing_self, names(nodup_indicies[i]))
  write(paste(names(nodup_indicies[i]), "NO MATCH", sep = "-match-"), file = 'testfiles.txt', append = T)
}

which(colnames(mixed_data_only) %in% "MSSM_RNA_PFC_356")
colnames(mixed_data_only)[248]
stripchart(mixed_data_only[,248], main = "MSSM_RNA_PFC_356")
which(mixed_data_only[,248] < 0.4)
which(rownames(mixed_data_only) %in% "MSSM_356")
rownames(mixed_data_only)[272]
stripchart(mixed_data_only[272,], main = "MSSM_356")
which(mixed_data_only[272,] < 0.4)

which(colnames(mixed_data_only) %in% "MSSM_RNA_PFC_154")
colnames(mixed_data_only)[73]
stripchart(mixed_data_only[,73], main = "MSSM_RNA_PFC_154")
which(mixed_data_only[,73] < 0.4)
which(rownames(mixed_data_only) %in% "MSSM_154")
rownames(mixed_data_only)[95]
stripchart(mixed_data_only[95,], main = "MSSM_154")
which(mixed_data_only[95,] < 0.4)

which(colnames(mixed_data_only) %in% "PENN_RNA_PFC_51")
colnames(mixed_data_only)[396]
stripchart(mixed_data_only[,396], main = "PENN_RNA_PFC_51")
which(mixed_data_only[,396] < 0.5)
which(rownames(mixed_data_only) %in% "PENN_51")
rownames(mixed_data_only)[427]
stripchart(mixed_data_only[427,], main = "PENN_51")
which(mixed_data_only[427,] < 0.5)

which(colnames(mixed_data_only) %in% "MSSM_RNA_PFC_21")
colnames(mixed_data_only)[121]
stripchart(mixed_data_only[,121], main = "MSSM_RNA_PFC_21")
which(mixed_data_only[,121] < 0.4)
mixed_data_only[c(239, 711), 121]
which(rownames(mixed_data_only) %in% "MSSM_BP_21")
rownames(mixed_data_only)[711]
stripchart(mixed_data_only[711,], main = "MSSM_BP_21")
which(mixed_data_only[711,] < 0.4)
