# KKD for Sage Bionetworks
# Jan. 14, 2014

library('synapseClient')
synapseLogin()
setwd('~/Computing/commonmind/analysis')
source('~/Computing/rnaseq_analysis_functions.R')

dataFile = synGet('syn2340130')
geneCounts  = read.delim(dataFile@filePath, header=T, row.names = 1)


##############
# Make final dataset, with bad samples excluded
##############

# which samples are repeated in this file
hist(table(metadata$DLPFC_RNA_isolation..Sample.RNA.ID))
which(table(metadata$DLPFC_RNA_isolation..Sample.RNA.ID) > 1)

# exclude the data marked for exclusion
colnames(metadata)
table(metadata$DLPFC_RNA_report..Exclude.)
table(metadata$DLPFC_RNA_report..Exclude.Reason)

metadata_freeze = metadata[-which(metadata$DLPFC_RNA_report..Exclude. == 1),]
table(metadata_freeze$DLPFC_RNA_report..Exclude.)
metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID = factor(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)
levels(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)
hist(table(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID))

# make a geneCounts matrix with reduced dataset
metadata_freeze = metadata_freeze[order(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID),]
geneCounts_freeze = geneCounts[,which(as.factor(colnames(geneCounts)) %in% metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)]
geneCounts_freeze = geneCounts_freeze[,order(as.factor(colnames(geneCounts_freeze)))]

head(colnames(geneCounts_freeze))
head(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)
tail(colnames(geneCounts_freeze))
tail(metadata_freeze$DLPFC_RNA_isolation..Sample.RNA.ID)



###############
# Detection analysis of biotypes
##############


head(rownames(geneCounts_freeze))
geneCounts_freeze_temp = geneCounts_freeze
geneCounts_freeze = addBiotype(geneCounts_freeze_temp)
rm(geneCounts_freeze_temp)

countDetectedByBiotype=function(sampleCounts,biotypes,filter=2){
  tofilter = which(sampleCounts < filter)
  tofilter = c(tofilter, which(sampleCounts == NA))
  filteredCounts = as.matrix(table(biotypes[-tofilter]))
  bioCounts = matrix(0, nrow = length(table(biotypes)), ncol = 1)
  rownames(bioCounts) = rownames(table(biotypes))
  bioCounts[which(rownames(bioCounts) %in% rownames(filteredCounts)),1] = filteredCounts[na.omit(match(rownames(bioCounts), rownames(filteredCounts))),1]
  return(bioCounts[,1])
}

totalBiotypeCounts = table(getBM(attributes=c("ensembl_gene_id", "gene_biotype"), mart=Hs)[,2])




plotDetectionByBiotype=function(inBiotypeCounts, inNormCounts){
  par(mfrow = c(3,3))
  corrStats = matrix(0, nrow = nrow(inNormCounts), ncol = 2)
  for (i in 1:27) {
     results = cor.test(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads), inBiotypeCounts[i,])
     corrStats[i,1] = as.vector(results[[3]])
     corrStats[i,2] = as.vector(results[[4]])
     
    plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads), inNormCounts[i,], main = rownames(inBiotypeCounts)[i], xlab = "mapped reads", ylab = "fraction detected", ylim = range(0,1))
    lines(lowess(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads), inNormCounts[i,]), col = "blue", lwd = 2)    
  }
  
  # interesting categories
  par(mfrow = c(2,3))
  interesting = c(2, 10,16,17,18,21)
  for (i in 1:length(interesting)) {
    plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads), inNormCounts[interesting[i],], main = paste(rownames(inBiotypeCounts)[interesting[i]], "\nn=", detectedTotal[interesting[i]], sep = ""), xlab = "mapped reads", ylab = "fraction detected", ylim = range(0,1))
    lines(lowess(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads), inNormCounts[interesting[i],]), col = "blue", lwd = 2)
  }
  par(mfrow = c(1,1))
  print(corrStats)
  return(corrStats)
}

regressDetectionOnMetadata=function(inDetectCounts, inMetaData){
  lmStats = matrix(0, nrow = nrow(inDetectCounts), ncol = 2)
  for (k in 1:nrow(inDetectCounts)) {
    results = summary(lm(inDetectCounts[k,] ~ inMetaData))
    lmStats[k,1] = as.vector(results$coefficients[2,1])
    lmStats[k,2] = as.vector(results$coefficients[2,4])
  }
  return(lmStats)
}

# detection - varying thresholds
corrResults = matrix(0, nrow = length(detectedTotal), ncol = 8)
lmResults = matrix(0, nrow = length(detectedTotal), ncol = 8)
j = 1
for (i in 1:4) {
  biotypeCountsPerSampleFilter = apply(geneCounts_freeze[,1:620], 2, countDetectedByBiotype, biotypes = geneCounts_freeze$biotype, filter=i)
  detectedTotal = totalBiotypeCounts[which(names(totalBiotypeCounts) %in% rownames(biotypeCountsPerSampleFilter))]
  normBiotypeCounts = apply(biotypeCountsPerSampleFilter, MARGIN=2, function(x) x/detectedTotal)
  corrResults[,j:(j+1)] = plotDetectionByBiotype(inBiotypeCounts=biotypeCountsPerSampleFilter,inNormCounts=normBiotypeCounts)
  lmResults[,j:(j+1)] = regressDetectionOnMetadata(inDetectCounts=normBiotypeCounts, inMetaData=log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads))
  j = j + 2
}


###############
# Detection analysis of expression levels
##############
## need to re-do this with length normalization so that expression level reflects expression only, not length
gene_lengths = getGeneLengths()
hist(log10(gene_lengths$length))


library('edgeR')
table(metadata_freeze$Dx)
#data.dge = DGEList(counts=geneCounts_freeze[,1:620],group=factor(metadata_freeze$Dx), remove.zeros=TRUE)
#data.dge = calcNormFactors(data.dge)
#save(data.dge,file="genecounts_freeze_phaseI.bz2", compress="bzip2")
load("genecounts_freeze_phaseI.bz2")
hist(cpm(data.dge,log=TRUE,normalized.lib.sizes=TRUE), main = "all samples, all genes", col = "thistle")
hist(aveLogCPM(data.dge, normalized.lib.sizes=TRUE), main = "gene averages across samples", col = "thistle", breaks = 30)

# control samples only with CPM
control.dge = data.dge[,which(data.dge$samples$group == "Control")]
hist(cpm(control.dge,log=TRUE,normalized.lib.sizes=TRUE), main = "control samples, all genes", col = "thistle")
control_gene_means = aveLogCPM(control.dge, normalized.lib.sizes=TRUE)
hist(control_gene_means, main = "gene averages across samples\ncontrol samples only", col = "thistle", breaks = 30, xlim = range(-5,15))


## control samples only with RPKM
hist(rpkm(control.dge, gene.length=gene_lengths$length[match(rownames(getCounts(control.dge)), gene_lengths$ensembl_gene_id)],log=TRUE,normalized.lib.sizes=TRUE), main = "control samples, all genes", col = "thistle", xlab = "RPKM")
control_gene_means = rowMeans(rpkm(control.dge, gene.length=gene_lengths$length[match(rownames(getCounts(control.dge)), gene_lengths$ensembl_gene_id)],log=TRUE,normalized.lib.sizes=TRUE))
hist(control_gene_means, main = "gene averages across samples\ncontrol samples only", col = "thistle", breaks = 50)

# separate datasets into high/med/low expression
range(control_gene_means, na.rm = TRUE)
lowExpress = which(control_gene_means < -7)
highExpress = which(control_gene_means > -2)
medExpress_temp = which(control_gene_means >= -7)
medExpress = setdiff(medExpress_temp, highExpress)
sum(length(lowExpress), length(highExpress), length(medExpress))
rm(medExpress_temp)

control_high.dge = control.dge[highExpress,]
control_med.dge = control.dge[medExpress,]
control_low.dge = control.dge[lowExpress,]


countDetected=function(x, filter=2){ length(which(x > filter)) }

control_high_detect = apply(getCounts(control_high.dge), MARGIN=2, FUN=countDetected)
control_med_detect = apply(getCounts(control_med.dge), MARGIN=2, FUN=countDetected)
control_low_detect = apply(getCounts(control_low.dge), MARGIN=2, FUN=countDetected)


# plot relationships - low expressors
#plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_low_detect, main = "detection of low-expression genes", xlab = "log10(mapped reads)", ylab = "number detected")
plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_low_detect/length(lowExpress), main = paste("detection of low-expression genes", "\nn=", length(lowExpress), sep = ""), xlab = "log10(mapped reads)", ylab = "fraction detected", ylim = range(0,1))
lines(lowess(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_low_detect/length(lowExpress)), col = "blue", lwd = 2)    
cor.test(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_low_detect)

# medium
#plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_med_detect, main = "detection of med-expression genes", xlab = "log10(mapped reads)", ylab = "number detected")
plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_med_detect/length(medExpress), main = paste("detection of med-expression genes", "\nn=", length(medExpress), sep = ""), xlab = "log10(mapped reads)", ylab = "fraction detected", ylim = range(0,1))
lines(lowess(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_med_detect/length(medExpress)), col = "blue", lwd = 2)    
cor.test(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_med_detect)


# high
#plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_high_detect, main = "detection of high-expression genes", xlab = "log10(mapped reads)", ylab = "number detected")
plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_high_detect/length(highExpress), main = paste("detection of high-expression genes", "\nn=", length(highExpress), sep = ""), xlab = "log10(mapped reads)", ylab = "fraction detected", ylim = range(0,1))
lines(lowess(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_high_detect/length(highExpress)), col = "blue", lwd = 2)    
cor.test(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)[which(metadata_freeze$Dx == "Control")], control_high_detect)

#############
# What genes are expressed at hi/med/low levels? Do they look like neurons?
##############

## trying to look at categories in each group
pcGenes = getByBiotype()
highExpress_ENSGpc = intersect(rownames(getCounts(control_high.dge)), pcGenes[,1])
medExpress_ENSGpc = intersect(rownames(getCounts(control_med.dge)), pcGenes[,1])
lowExpress_ENSGpc = intersect(rownames(getCounts(control_low.dge)), pcGenes[,1])
write.csv(lowExpress_ENSGpc, file = "lowexpress_ENSGpc", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(highExpress_ENSGpc, file = "highexpress_ENSGpc", quote = FALSE, row.names = FALSE, col.names = FALSE)

# probably need to attempt by taking known groups of genes and looking for what fraction of them fall in different groups.

printFractionalCounts=function(categoryENSGvec, categoryName){
  print(paste("Number of genes in category", length(categoryENSGvec), sep = " "))
  print(paste("High expression genes as fraction of all genes:", length(highExpress) / nrow(control.dge), sep = " "))
  print(paste("Fraction of", categoryName, "genes found in high expression group:", length(intersect(rownames(getCounts(control_high.dge)), categoryENSGvec))/length(categoryENSGvec), sep = " "))
  
  print(paste("Medium expression genes as fraction of all genes:", length(medExpress) / nrow(control.dge), sep = " "))
  print(paste("Fraction of", categoryName, "genes found in med expression group:", length(intersect(rownames(getCounts(control_med.dge)), categoryENSGvec))/length(categoryENSGvec), sep = " "))
  
  print(paste("Low expression genes as fraction of all genes:", length(lowExpress) / nrow(control.dge), sep = " "))
  print(paste("Fraction of", categoryName, "genes found in low expression group:", length(intersect(rownames(getCounts(control_low.dge)), categoryENSGvec))/length(categoryENSGvec), sep = " "))  
}

# genes in neuron-related cellular component GO
neuron_part = getGenesForGOOffspring("GO:0097458")
synapse = getGenesForGOOffspring("GO:0045202")
vesicle = getGenesForGOOffspring("GO:1990008") # nothing found

neuron_and_synape = union(neuron_part, synapse)
printFractionalCounts(neuron_and_synape, "neuron-synapse")

# genes in neuron birth-growth-death biological process GO
neurogenesis = getGenesForGOOffspring("GO:0022008",go="BP")	
neuron_death = getGenesForGOOffspring("GO:0070997", go="BP")
neuron_diff = getGenesForGOOffspring("GO:0030182", go="BP")
neuron_bgd = union(union(neurogenesis, neuron_death), neuron_diff)

printFractionalCounts(neuron_bgd, "neuron birth-growth-death")

length(intersect(neuron_and_synape, neuron_bgd))
# not sure what other categories to do, or if I should expand BP and MF more. CC seems to show well that we are getting neuron=related expression. I could contrast with another CC?



#######
# Detection by biotype and expression?
#######

# run lm on detection per biotype
# sets of two columns: first is beta estimate for slope, second is p-value
print(lmResults)

# run ANOVA/lm on detection per biotype and expression level
# sets of two columns: first is beta estimate for slope, second is p-value


##########
# How do the expression distributions of biotypes compare? Are lincRNAs hi/med/lo compared to protein coding?
# plot distn of expression levels per biotype, length-normalized.
##########


#####
# lncRNAs
####

lincRNA_geneCounts = geneCounts_freeze[which(geneCounts_freeze$biotype == "lincRNA"),]

# get length
hist(log10(gene_lengths[which(gene_lengths$ensembl_gene_id %in% rownames(lincRNA_geneCounts)),4]), main = "lincRNA gene lengths", xlab = "log10(gene end - gene start)", breaks = 20)

# what is range of length?
linc_lengths = gene_lengths[which(gene_lengths$ensembl_gene_id %in% rownames(lincRNA_geneCounts)),c(1,4)]
range(linc_lengths[,2])

# detection by length -- as we add seq depth, are we detecting more longer or shorter lincRNAs
llincRNAs = which(linc_lengths[,2] > 1e5)
slincRNAs = which(linc_lengths[,2] <= 1e5)

llincRNA_detect = apply(lincRNA_geneCounts[llincRNAs,1:620], MARGIN=2, FUN=countDetected)
slincRNA_detect = apply(lincRNA_geneCounts[slincRNAs,1:620], MARGIN=2, FUN=countDetected)

hist(llincRNA_detect)
hist(slincRNA_detect)

# note this includes all samples, not just control samples
# make another plot where fraction is out of total lincRNAs?
plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads), llincRNA_detect[1:620]/length(llincRNAs), main = paste("detection of long lincRNAs", "\nn=", length(llincRNAs), sep = ""), xlab = "log10(mapped reads)", ylab = "fraction detected out of total for this dataset", ylim = range(0,1))
summary(lm(llincRNA_detect/length(llincRNAs) ~ log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)))


plot(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads), slincRNA_detect/length(slincRNAs), main = paste("detection of short-med lincRNAs", "\nn=", length(slincRNAs), sep = ""), xlab = "log10(mapped reads)", ylab = "fraction detected out of total for this dataset", ylim = range(0,1))
lines(lowess(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads), slincRNA_detect/length(slincRNAs)), col = "blue", lwd = 2)    
cor.test(log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads),slincRNA_detect/length(slincRNAs))
summary(lm(slincRNA_detect/length(slincRNAs) ~ log10(metadata_freeze$DLPFC_RNA_report..Mapped.Reads)))

