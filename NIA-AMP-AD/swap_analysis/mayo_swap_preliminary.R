library('synapseClient')
library('edgeR')
synapseLogin()

# Get gene counts data and metadata
mayo = synGet('syn5553124')
mayoData = read.table(getFileLocation(mayo),row.names = 1)
mssm = synGet('syn5553119')
mssmData = read.table(getFileLocation(mssm), row.names = 1)


metadataEnt = synGet("syn3817650")
metadata = read.csv(getFileLocation(metadataEnt))

# Fix the one mayo filename that doesn't have "TCX" in name
colnames(mssmData)
colnames(mayoData)
z = which(colnames(mayoData) == "X11427")
colnames(mayoData)[z] = "X11427_TCX"

# Truncate the gene names at the "."
x = sapply(strsplit(rownames(mssmData), "[.]"), function(x){ return(x[1]) } )
rownames(mssmData) = x
x = sapply(strsplit(rownames(mayoData), "[.]"), function(x){ return(x[1]) } )
rownames(mayoData) = x


# Normalize for read depth
mayoR.dge = DGEList(counts=mayoData,remove.zeros=TRUE)
mayoR.dge = calcNormFactors(mayoR.dge)
mayoR_counts = cpm(mayoR.dge,log = TRUE,normalized.lib.sizes = TRUE)

mssmR.dge = DGEList(counts=mssmData,remove.zeros = TRUE)
mssmR.dge = calcNormFactors(mssmR.dge)
mssmR_counts = cpm(mssmR.dge,log = TRUE,normalized.lib.sizes = TRUE)

mayoPalx.dge = calcNormFactors(filterByFractionPresent(inCounts = mayoR.dge,fraction = 0.6))
mssmPalx.dge = calcNormFactors(filterByFractionPresent(inCounts = mssmR.dge,fraction = 0.6))

mayoCountsPalx = cpm(mayoPalx.dge,log = TRUE,normalized.lib.sizes = TRUE)
mssmCountsPalx = cpm(mssmPalx.dge,log = TRUE,normalized.lib.sizes = TRUE)

