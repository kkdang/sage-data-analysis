library('synapseClient')
library('edgeR')
synapseLogin()

# Get gene counts data and metadata
mayo = synGet('syn5583761')
mayoData = read.table(getFileLocation(mayo),row.names = 1)
mssm = synGet('syn5583762')
mssmData = read.table(getFileLocation(mssm), row.names = 1)

# Truncate the gene names at the "."
x = sapply(strsplit(rownames(mssmData), "[.]"), function(x){ return(x[1]) } )
rownames(mssmData) = x
x = sapply(strsplit(rownames(mayoData), "[.]"), function(x){ return(x[1]) } )
rownames(mayoData) = x


metadataEnt = synGet("syn3205337")
metadata = read.delim(getFileLocation(metadataEnt))

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

rm(mssmData, mayoData)

mayoPalx.dge = calcNormFactors(filterByFractionPresent(inCounts = mayoR.dge,fraction = 0.6))
mssmPalx.dge = calcNormFactors(filterByFractionPresent(inCounts = mssmR.dge,fraction = 0.6))

mayoCountsPalx = cpm(mayoPalx.dge,log = TRUE,normalized.lib.sizes = TRUE)
mssmCountsPalx = cpm(mssmPalx.dge,log = TRUE,normalized.lib.sizes = TRUE)
