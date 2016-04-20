## QC AMP-AD re-processed count files
## KKD for Sage Bionetworks
## April 5, 2016



library('rGithubClient')
source('/Users/kristen/Computing/external_software/rgithubclient_authenticate.R')
library(edgeR)
library(gplots)
library(RColorBrewer)
library(synapseClient)
sourceRepoFile(sageCode, "rnaseq_analysis_functions.R")
synapseLogin()

# featurecounts assignment metrics
plotFCmetrics=function(synID){
  fcAssignMetrics = read.delim(getFileLocation(synGet(synID)))
  head(fcAssignMetrics)
  fc.m = data.matrix(fcAssignMetrics)
  head(fc.m)
  fcTotals = colSums(fc.m)
  
  op = par(mfrow = c(3,3))
  for (i in 2:nrow(fc.m)){
    hist(as.numeric(fc.m[i,]), col = "honeydew", main = rownames(fc.m)[i], xlab = "counts")
  }
  
  for (i in 2:nrow(fc.m)){
    hist(as.numeric(fc.m[i,])/fcTotals, col = "honeydew", main = rownames(fc.m)[i], xlab = "fractional allocation")
  }
  par(op)
}

countFileList = synQuery('select id,name from file where parentId=="syn5870306"')
countFileList
summaryFiles = countFileList$file.id[grep("summary", countFileList$file.name)]
sapply(summaryFiles, plotFCmetrics)


ml = read.delim(getFileLocation(synGet("syn5870315")))
mt = read.delim(getFileLocation(synGet("syn5870316")))
sl = read.delim(getFileLocation(synGet("syn5870313")))
st = read.delim(getFileLocation(synGet("syn5870314")))
rl = read.delim(getFileLocation(synGet("syn5870317")))
rt = read.delim(getFileLocation(synGet("syn5870318")))

hist(as.numeric(ml[6,]) - as.numeric(mt[6,]))
hist(as.numeric(sl[6,]) - as.numeric(st[6,]))
range(as.numeric(sl[6,]) - as.numeric(st[6,]))
hist(as.numeric(rl[6,]) - as.numeric(rt[6,]))

which(rl[6,] == rt[6,1])
which(as.numeric(rl[6,]) == as.numeric(rt[6,1]))

#  biotypes profile and ribosomal fraction

biotypeAndRibosomal=function(synID){
  dataEnt = read.delim(getFileLocation(synGet(synID)))
  rownames(dataEnt) = b = sapply(rownames(dataEnt),function(x){unlist(strsplit(as.character(x), split = "[.]"))[1]})
  data_wbiotypes = addBiotype(inCounts = dataEnt)
  pie(table(data_wbiotypes$biotype)) # detected biotypes, not fractional allocation
  data_ribosomal = data_wbiotypes[which(data_wbiotypes$biotype == "rRNA"),1:9]
  hist(colSums(data_ribosomal) / colSums(data_wbiotypes[,1:9]))
  return(data_wbiotypes)
}

countFiles = countFileList$file.id[grep("summary", countFileList$file.name,invert = TRUE)]
sapply(countFiles, biotypeAndRibosomal)


# correlation heatmaps

op = par(mar = c(12,4,4,8))
for (i in 1:length(countFiles)) {
  synID = countFiles[i]
  dataEnt = read.delim(getFileLocation(synGet(synID)))
  data.dge = DGEList(counts=dataEnt,remove.zeros=TRUE)
  data.dge = calcNormFactors(data.dge)
  data_counts = cpm(data.dge,log = TRUE,normalized.lib.sizes = TRUE)
  data.dend = heatmap.2(cor(data_counts,method = "spearman"),dendrogram = "col",scale = "none",main = synID,trace = "none",breaks = seq(0.80,1,0.001), margins = c(7,7), labRow = NULL, labCol = NULL)
}
