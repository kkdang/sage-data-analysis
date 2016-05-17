##### Mt. Sinai data ######## 
# RNAseq

library(synapseClient)
synapseLogin()

commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="geneExp", fileType="bam")

datasetAnnotations = list(study="MSBB",center="MSBB",modelSystem=FALSE,organism="HomoSapiens",platform="IlluminaHiSeq2500")

# Batch 1
bamlist = synQuery("select id, name from file where parentId=='syn4055270'")
onlyBam = bamlist[grep(pattern = "bam",bamlist$file.name),]
onlyFq = bamlist[grep(pattern = "unmapped.fq",bamlist$file.name),]
for (i in 1:nrow(onlyBam)){
  temp = synGet(onlyBam$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}
commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="geneExp", fileType="fastq")
for (i in 1:nrow(onlyFq)){
  temp = synGet(onlyFq$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}

# Batch 2
commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="geneExp", fileType="bam")
bamlist = synQuery("select id, name from file where parentId=='syn4883087'")
onlyBam = bamlist[grep(pattern = "bam",bamlist$file.name),]
onlyFq = bamlist[grep(pattern = "unmapped.fq",bamlist$file.name),]
for (i in 1:nrow(onlyBam)){
  temp = synGet(onlyBam$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}
commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="geneExp", fileType="fastq")
for (i in 1:nrow(onlyFq)){
  temp = synGet(onlyFq$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}

# Batch 3
commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="geneExp", fileType="bam")
bamlist = synQuery("select id, name from file where parentId=='syn5480006'")
onlyBam = bamlist[grep(pattern = "bam",bamlist$file.name),]
onlyFq = bamlist[grep(pattern = "unmapped.fq",bamlist$file.name),]
for (i in 1:nrow(onlyBam)){
  temp = synGet(onlyBam$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}
commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="geneExp", fileType="fastq")
for (i in 1:nrow(onlyFq)){
  temp = synGet(onlyFq$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}

