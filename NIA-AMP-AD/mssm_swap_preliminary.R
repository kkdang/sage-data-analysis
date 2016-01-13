library('synapseClient')
library('edgeR')
synapseLogin()

# Get gene counts data and metadata
mayo = synGet('syn5583761')
mayoData = read.table(getFileLocation(mayo),row.names = 1)
mssm = synGet('syn5583762')
mssmData = read.table(getFileLocation(mssm), row.names = 1)


metadataEnt = synGet("syn3205337")
metadata = read.delim(getFileLocation(metadataEnt))

# Fix the one mayo filename that doesn't have "TCX" in name
colnames(mssmData)
match(colnames(mssmData), metadata$Sample.ID)
colnames(mayoData)
#x = paste("BM_10_", substr(colnames(mayoData),start = 4,6), sep = "")
#colnames(mayoData) = x
colnames(mayoData)
colnames(mssmData)


# Normalize for read depth
mayoR.dge = DGEList(counts=mayoData,remove.zeros=TRUE)
mayoR.dge = calcNormFactors(mayoR.dge)
mayoR_counts = cpm(mayoR.dge,log = TRUE,normalized.lib.sizes = TRUE)

mssmR.dge = DGEList(counts=mssmData,remove.zeros = TRUE)
mssmR.dge = calcNormFactors(mssmR.dge)
mssmR_counts = cpm(mssmR.dge,log = TRUE,normalized.lib.sizes = TRUE)
