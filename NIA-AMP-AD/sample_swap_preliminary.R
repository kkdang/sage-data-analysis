library('synapseClient')
library('edgeR')
synapseLogin()

# Get gene counts data and metadata
mayo = synGet('syn4904891')
mayoData = read.csv(getFileLocation(mayo),row.names = 1)
mssm = synGet('syn4904889')
mssmData = read.csv(getFileLocation(mssm), row.names = 1)
broad = synGet('syn4904859')
broadData = read.csv(getFileLocation(broad), row.names = 1)


metadataEnt = synGet('syn4166661')
metadata = read.csv(getFileLocation(metadataEnt),header = TRUE)
metadataEnt2 = synGet('syn4300313')
metadata2 = read.delim(getFileLocation(metadataEnt2))

# Change IDs in Broad data to match the other 2.
# IDs in MSSM and Mayo data files and metadata are ProjID
colnames(mssmData)
colnames(mayoData)

# IDs from ROSMAP are SampleID
head(colnames(broadData))
mrna_id_strip = sapply(as.character(metadata2$ID),function(z){paste(unlist(strsplit(z,split = "_"))[1:2],collapse = "_")})
colnames(broadData) = paste("X",metadata2$projid[match(colnames(broadData), paste("X",mrna_id_strip,sep = ''))],sep = "")
colnames(broadData)


# Normalize for read depth
mayoR.dge = DGEList(counts=mayoData,remove.zeros=TRUE)
mayoR.dge = calcNormFactors(mayoR.dge)
mayoR_counts = cpm(mayoR.dge,log = TRUE,normalized.lib.sizes = TRUE)

mssmR.dge = DGEList(counts=mssmData,remove.zeros = TRUE)
mssmR.dge = calcNormFactors(mssmR.dge)
mssmR_counts = cpm(mssmR.dge,log = TRUE,normalized.lib.sizes = TRUE)

broadR.dge = DGEList(counts=broadData,remove.zeros = TRUE)
broadR.dge = calcNormFactors(broadR.dge)
broadR_counts = cpm(broadR.dge,log = TRUE,normalized.lib.sizes = TRUE)