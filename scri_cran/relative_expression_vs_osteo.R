#! /usr/bin/env Rscript
# Jan. 30, 2017
# KKD for Sage Bionetworks

library(synapseClient)
synapseLogin()

library('githubr')
source('/Users/kkdang/Computing/rgithubclient_authenticate.R')
sourceRepoFile(sageCode, "scri_cran/process_metadata_validation.R")
setwd('~/Computing/cranio/')
library(biomaRt)
library(gplots)
library(edgeR)

Hs = useMart('ENSEMBL_MART_ENSEMBL')
Hs = useDataset(dataset = "hsapiens_gene_ensembl",mart = Hs)


# get residual data
resid = read.delim(getFileLocation(synGet("syn5873539")))

# get transcript lengths
transLengths = getBM(attributes = c("ensembl_gene_id", "transcript_length", "ensembl_transcript_id"),filters = "transcript_biotype",values = "protein_coding",mart = Hs)
head(transLengths)

getLongestTranscript=function(ensg,transLengthsBM){
  geneLines = which(transLengthsBM$ensembl_gene_id %in% ensg)
  return(max(transLengthsBM$transcript_length[geneLines]))
}
longLengths = lapply(as.list(unique(transLengths$ensembl_gene_id)),function(s){getLongestTranscript(ensg = s,transLengths)})
names(longLengths) = unique(transLengths$ensembl_gene_id)

# get gene lists
osteoAll = read.delim('~/Downloads/GeneList-from QiagenRT2-profiler_for Kristin to check expression in candidates.txt',header = FALSE)
osteo = data.frame(geneNames = osteoAll[-which(duplicated(osteoAll[,1])),])
symb = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),filters = "hgnc_symbol",values = osteo[,1],mart = Hs)
osteo$ENSG = symb$ensembl_gene_id[match(osteo[,1], symb$hgnc_symbol)]
head(osteo)

ofInterest = data.frame(geneNames = c("AXL","PDGFRA", "PIEZO1", "FLNA", "FLNB", "FLNC"))
symb = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),filters = "hgnc_symbol",values = ofInterest$geneNames,mart = Hs)
ofInterest$ENSG = symb$ensembl_gene_id[match(ofInterest$geneNames, symb$hgnc_symbol)]
head(ofInterest)

# calculate median transcript length and get correction factor for each gene
ofInterest$transLength = longLengths[match(as.character(ofInterest$ENSG), names(longLengths))]
osteo$transLength = longLengths[match(as.character(osteo$ENSG), names(longLengths))]

medOfAll = median(c(as.numeric(ofInterest$transLength), as.numeric(osteo$transLength)))

ofInterest$adj = medOfAll/as.numeric(ofInterest$transLength)
osteo$adj = medOfAll/as.numeric(osteo$transLength)

head(osteo)
head(ofInterest)

# divide genes by correction factor 
ofInterestExp = resid[which(rownames(resid) %in% ofInterest$ENSG),]
ofInterest = ofInterest[match(rownames(ofInterestExp),ofInterest$ENSG),]
rownames(ofInterestExp)
ofInterest$ENSG
ofInterestExp = ofInterestExp * ofInterest$adj


osteoExp = resid[which(rownames(resid) %in% osteo$ENSG),]
toRemove = setdiff(osteo$ENSG, rownames(osteoExp))
osteoTrimmed = osteo[-which(osteo$ENSG %in% toRemove),]
rm(osteo)
osteoTrimmed = osteoTrimmed[match(rownames(osteoExp),osteoTrimmed$ENSG),]
head(osteoTrimmed$ENSG)
head(rownames(osteoExp))

osteoExp = osteoExp * osteoTrimmed$adj
head(osteoExp)

# heatmap non-scaled of geneList vs osteo genes across all patients
combined = rbind(osteoExp,ofInterestExp)
heatmap.2(data.matrix(combined),scale = NULL,col = bluered(9),trace = "none",RowSideColors = c(rep("white",nrow(osteoExp)),rep("orange", nrow(ofInterestExp))),labRow = "")

# boxplot of corrected gene values for gene list and osteo genes.
boxplot(t(combined), las = 2)



## Another look, using CPM data instead.
counts = read.csv(getFileLocation(synGet("syn2820309")),row.names = 1)
head(counts)
counts.dge =DGEList(counts = counts,remove.zeros = TRUE) 
counts.dge = calcNormFactors((counts.dge))


shortCounts.dge = counts.dge[which(rownames(counts.dge) %in% rownames(combined)),]
shortCounts.dge = shortCounts.dge[match(rownames(combined),rownames(shortCounts.dge)),]
head(rownames(shortCounts.dge))
head(rownames(combined))
head(osteoTrimmed$ENSG)
counts.rpkm = rpkm(shortCounts.dge,gene.length = as.numeric(as.vector(c(osteoTrimmed$transLength,ofInterest$transLength))),normalized.lib.sizes = TRUE,log = TRUE)


heatmap.2(counts.rpkm,scale = NULL,col = bluered(9),trace = "none",RowSideColors = c(rep("white",nrow(osteoExp)),rep("orange", nrow(ofInterestExp))),labRow = "")
boxplot(counts.rpkm, las = 2,use.cols = FALSE,col = c(rep("white",nrow(osteoExp)),rep("orange", nrow(ofInterestExp))),names = c(as.character(osteoTrimmed$geneNames), as.character(ofInterest$geneNames)),cex.axis = 0.8, ylab = "RPKM")
