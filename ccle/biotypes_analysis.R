library("biomaRt")

setwd('~/Computing//H3-ccle')
data = read.delim('htseq-ccle.log', sep = "\t", header = F)
totalcounts = sum(data[,2])

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

attributes = listAttributes(ensembl)
y = grep('biotype', attributes[,1], ignore.case=T)
attributes[y,1]

filters = listFilters(ensembl)
x = grep('biotype', filters[,1], ignore.case=T)
filters[x,1]


biotypes = getBM(c('ensembl_gene_id','gene_biotype'), filters = c('ensembl_gene_id'), values=data[,1], mart=ensembl)
data[,3] = biotypes[match(data[,1], biotypes[,1]),2]

# distribution of detected types
pie(table(data[which(data[,2] > 0),3]))

# distribution of counts per detected type
types = names(table(data[,3]))
summed_types = rep(NA, length(types))
names(summed_types) = types
for (i in 1:length(types)) {
  summed_types[i] = sum(data[which(data[,3] == types[i]),2])
}
pie(summed_types)