# Dec. 13, 2015
# KKD for Sage Bionetworks
# Annotates BAM data for AMP-AD


library(synapseClient)
synapseLogin()

commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="geneExp", fileType="bam")



##### Mayo data ######## TCX #
datasetAnnotations = list(study="MayoRNAseq",center="UFL-Mayo-ISB",modelSystem=FALSE,tissueTypeAbrv="TCX",tissueType="temporalCortex",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn4894912'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}



##### Mayo data ######## Cerebellum #
datasetAnnotations = list(study="MayoRNAseq",center="UFL-Mayo-ISB",modelSystem=FALSE,tissueTypeAbrv="CBE",tissueType="cerebellum",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn5049322'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}

##### Mayo data ######## Sample swap seq at Sinai #
datasetAnnotations = list(study="sampleSwap",center="MSSM",modelSystem=FALSE,tissueTypeAbrv="TCX",tissueType="temporalCortex",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn3537578'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}


##### Mayo data ######## pilot PSP data #
datasetAnnotations = list(study="MayoPilot",center="UFL-Mayo-ISB",modelSystem=FALSE,tissueTypeAbrv="TCX",tissueType="temporalCortex",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn5584594'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}

##### Mayo data ######## pilot AD data #
datasetAnnotations = list(study="MayoPilot",center="UFL-Mayo-ISB",modelSystem=FALSE,tissueTypeAbrv="TCX",tissueType="temporalCortex",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn5580964'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}



##### ROSMAP data ######## 
# RNAseq
datasetAnnotations = list(study="ROSMAP",center="Broad-Rush",modelSystem=FALSE,tissueTypeAbrv="PFC",tissueType="dorsolateralPrefrontalCortex",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn4164376'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}

# ChIPseq
notSoCommonAnnotations = list(consortium="AMP-AD", dataType="DNA", assay="ChIPseq", fileType="bam")
datasetAnnotations = list(study="ROSMAP",center="Broad-Rush",modelSystem=FALSE,tissueTypeAbrv="PFC",tissueType="dorsolateralPrefrontalCortex",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn4896408'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(notSoCommonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}



##### Broad MDMi data ######## 
datasetAnnotations = list(study="BroadMDMi",center="Broad-Rush",modelSystem=TRUE,organism="HomoSapiens",tissueType="peripheralBloodMononuclearCell")

bamlist = synQuery("select id, name from file where parentId=='syn4228560'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}


##### Broad iPSC data ######## 
datasetAnnotations = list(study="BroadiPSC",center="Broad-Rush",modelSystem=TRUE,organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn4228582'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}


##### ROSMAP sample swap data ######## 
datasetAnnotations = list(study="sampleSwap",center="MSSM",modelSystem=FALSE,tissueTypeAbrv="PFC",tissueType="dorsolateralPrefrontalCortex",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn4164988'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}


datasetAnnotations = list(study="sampleSwap",center="UFL-Mayo-ISB",modelSystem=FALSE,tissueTypeAbrv="PFC",tissueType="dorsolateralPrefrontalCortex",organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn4904772'")
for (i in 1:nrow(bamlist)){
  temp = synGet(bamlist$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}



##### Mt. Sinai data ######## 
# RNAseq
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

##### MSSM iPSC data ######## 
commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="geneExp", fileType="bam")

datasetAnnotations = list(study="MSSMiPSC",center="MSSM",modelSystem=TRUE,organism="HomoSapiens")

bamlist = synQuery("select id, name from file where parentId=='syn5986894'")
onlyBam = bamlist[grep(pattern = "bam",bamlist$file.name),]
onlyFq = bamlist[grep(pattern = "unmapped",bamlist$file.name),]
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

mssmipsc = read.csv(getFileLocation(synGet("syn6040297")))
head(mssmipsc)
head(bamlist)
table(mssmipsc$Celltype)
n = which(mssmipsc$Celltype == "Neurons")
npsc = which(mssmipsc$Celltype == "NPCs")

sampleIds = sapply(as.list(bamlist$file.name), function(x) { unlist(strsplit(x, split = "[.]"))[1] })
bamsToEdit = which(sampleIds %in% mssmipsc$Library[n])
for (i in 1:length(bamsToEdit)){
  tempEnt = synGet(bamlist$file.id[bamsToEdit[i]], downloadFile = FALSE)
  synSetAnnotation(tempEnt, "celltype") = "iPSCderivedMatureNeuron"
  synStore(tempEnt, forceVersion = FALSE)
}

bamsToEdit = which(sampleIds %in% mssmipsc$Library[npsc])
for (i in 1:length(bamsToEdit)){
  tempEnt = synGet(bamlist$file.id[bamsToEdit[i]], downloadFile = FALSE)
  synSetAnnotation(tempEnt, "celltype") = "iPSCderivedNeuralProgenitorCell"
  synStore(tempEnt, forceVersion = FALSE)
}

##### Mayo data ######## TCX #
commonAnnotations = list(consortium="AMP-AD", dataType="mRNA", assay="RNAseq", dataSubType="differentialExpression", fileType="tsv")

datasetAnnotations = list(study="MayoRNAseq",center="UFL-Mayo-ISB",modelSystem=FALSE,tissueTypeAbrv="TCX",tissueType="temporalCortex",organism="HomoSapiens")


entityList = synQuery("select id, name from file where parentId=='syn6090803'")
for (i in 1:nrow(entityList)){
  temp = synGet(entityList$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}
entityList = synQuery("select id, name from file where parentId=='syn6090804'")
for (i in 1:nrow(entityList)){
  temp = synGet(entityList$file.id[i],downloadFile = FALSE)
  synSetAnnotations(temp) = c(commonAnnotations, datasetAnnotations)
  temp = synStore(temp,forceVersion = FALSE)
}
