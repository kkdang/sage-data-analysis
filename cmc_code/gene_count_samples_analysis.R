# KKD for Sage Bionetworks
# Jan. 14, 2014

library('synapseClient')
library('edgeR')
synapseLogin()
setwd('~/Computing/commonmind/analysis')

dataFile = synGet('syn2340130')
geneCounts  = read.delim(dataFile@filePath, header=T, row.names = 1)

##########
# Alignment metrics analysis
##########
metadata = read.csv(getFileLocation(syn_report))

# Alignment overview
par(mfrow = c(1,3))
hist(metadata$DLPFC_RNA_report..Mapped.Reads, col = "deeppink4")
hist(metadata$DLPFC_RNA_report..Genes.Detected, col = "deeppink4")
hist(metadata$DLPFC_RNA_report..Percent.Aligned, col = "deeppink")

# RIN effect
plot(metadata$DLPFC_RNA_isolation..RIN, log10(metadata$DLPFC_RNA_report..Mapped.Reads), main = "Alignment yield (log10) vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Mapped.Reads, main = "Alignment yield vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Percent.Aligned, main = "Fraction reads aligned vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, log10(metadata$DLPFC_RNA_report..rRNA.Rate), main = "rRNA aligned vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Intragenic.Rate, main = "Intragenic aligned vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Exonic.Rate, main = "Exonic aligned vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Intronic.Rate, main = "Intronic aligned vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Intergenic.Rate, main = "Intergenic aligned vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Genes.Detected, main = "Genes detected aligned vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Transcripts.Detected, main = "Transcripts detected aligned vs RIN")
plot(metadata$DLPFC_RNA_isolation..RIN, metadata$DLPFC_RNA_report..Expression.Profiling.Efficiency, main = "Exp profiling efficiency vs RIN")

# which samples are repeated in this file
hist(table(metadata$DLPFC_RNA_isolation..Sample.RNA.ID))
which(table(metadata$DLPFC_RNA_isolation..Sample.RNA.ID) > 1)

# exclude the data marked for exclusion
colnames(metadata)
table(metadata$DLPFC_RNA_report..Exclude.)
table(metadata$DLPFC_RNA_report..Exclude.Reason)

metadata_reduced = metadata[-which(metadata$DLPFC_RNA_report..Exclude. == 1),]
levels(metadata_reduced$DLPFC_RNA_isolation..Sample.RNA.ID)
metadata_reduced$DLPFC_RNA_isolation..Sample.RNA.ID = factor(metadata_reduced$DLPFC_RNA_isolation..Sample.RNA.ID)
levels(metadata_reduced$DLPFC_RNA_isolation..Sample.RNA.ID)
hist(table(metadata_reduced$DLPFC_RNA_isolation..Sample.RNA.ID))


##############
# sample-level variation
##############

table(metadata$Dx)

metadata_seq = metadata_reduced[which(metadata_reduced$DLPFC_RNA_isolation..Sample.RNA.ID %in% as.factor(colnames(geneCounts))),]
metadata_seq = metadata_seq[order(metadata_seq$DLPFC_RNA_isolation..Sample.RNA.ID),]

geneCounts = geneCounts[,order(as.factor(colnames(geneCounts)))]
head(colnames(geneCounts))
head(metadata_seq$DLPFC_RNA_isolation..Sample.RNA.ID)

tail(colnames(geneCounts))
tail(metadata_seq$DLPFC_RNA_isolation..Sample.RNA.ID)

data.dge = DGEList(counts=geneCounts,group=factor(metadata_seq$Dx))
data.dge = calcNormFactors(data.dge)



## MDS plots
bp.dge = data.dge[,which(data.dge$samples$group == "BP")]
scz.dge = data.dge[,which(data.dge$samples$group == "SCZ")]
control.dge = data.dge[,which(data.dge$samples$group == "Control")]

# Control samples
pdf(file="MDS_bcv_control_default.pdf")
plotMDS.DGEList(x=control.dge, method = "bcv", main = "control samples")
dev.off()

pdf(file="MDS_logFC_control_default.pdf")
plotMDS.DGEList(x=control.dge, main = "control samples")
dev.off()

# BP
pdf(file="MDS_bcv_BP_default.pdf")
plotMDS.DGEList(x=bp.dge, method = "bcv", main = "BP samples")
dev.off()

pdf(file="MDS_logFC_BP_default.pdf")
plotMDS.DGEList(x=bp.dge, main = "BP samples")
dev.off()

# SCZ
pdf(file="MDS_bcv_SCZ_default.pdf")
plotMDS.DGEList(x=scz.dge, method = "bcv", main = "SCZ samples")
dev.off()

pdf(file="MDS_logFC_SCZ_default.pdf")
plotMDS.DGEList(x=scz.dge, main = "SCZ samples")
dev.off()

# all samples
plotLabels = as.numeric(metadata_reduced$Dx)
plotColors = c("red", "green", "blue")
plotMDS.DGEList(x=data.dge, labels=plotLabels, col = plotColors[plotLabels])

plotMDS.DGEList(x=data.dge, labels=plotLabels, col = plotColors[plotLabels], method = "bcv")

## need global look at source of variation (Menachem did this already for poster)


## heatmap of sample-sample correlation



## kernel density plots
densityPlot(bp.dge, mainLabel="BP")
densityPlot(scz.dge, mainLabel="SCZ")
densityPlot(control.dge, mainLabel="control")


## quick look at group level differences
# bp and control are more similar to each other than to scz
bp.ave = aveLogCPM(bp.dge)
scz.ave = aveLogCPM(scz.dge)
control.ave = aveLogCPM(control.dge)

plot(density(bp.ave), col = "red", lwd = 2, xlim = range(-5,2), main = "per-sample ave density plot")
lines(density(scz.ave), col = "blue", lwd = 2)
lines(density(control.ave), col = "green", lwd = 2)
legend("topright", legend=c("bp", "control", "schiz"), col=plotColors,pch=16)


y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)