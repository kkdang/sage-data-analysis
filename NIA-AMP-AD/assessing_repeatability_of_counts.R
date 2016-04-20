# Looking at repeatability of featurecounts output

library(synapseClient)
synapseLogin()

numDiffs = read.delim("~/Computing/NIA-AMP-AD/reprocessed/fc_diffs/lines_per_sample", header = FALSE, sep = " ")

# number of differences per sample
hist((numDiffs[,1]/4)-1, col = "honeydew", breaks = 40, xlab = "number of genes different", main = "genes counted differently\nper sample")

# fraction of samples having more than one gene different
print(sum(numDiffs[,1] > 8) / nrow(numDiffs))


# magnitude of differences
mag = read.delim(getFileLocation(synGet("syn5964656")), header = FALSE)
head(mag)
hist(mag$V3, col = "honeydew", xlab = "difference (in counts)", main = "size of differences between counts\nacross all differences")
hist(abs(mag$V3)/mag$V4, col = "honeydew", breaks = 20)


# per-sample mean magnitudes
allSampNames = levels(mag[,1])
hist(sapply(as.list(allSampNames),function(x) { mean(mag$V3[mag[,1] == x]) }), col = "honeydew", xlab = "mean difference across genes", main = "per-sample average differences")
hist(sapply(as.list(allSampNames),function(x) { mean(mag$V3[mag[,1] == x]) }), col = "honeydew", xlab = "mean difference across genes", main = "per-sample average differences", xlim = range(-500,5000), breaks = 200)
hist(sapply(as.list(allSampNames),function(x) { mean(mag$V3[mag[,1] == x]) }), col = "honeydew", xlab = "mean difference across genes", main = "per-sample average differences", xlim = range(-50,100), breaks = 2000)


trunSampNames = sapply(as.list(allSampNames),function(x){ substr(x, 1, 10) })

boxplot(sapply(as.list(allSampNames),function(x) { mag$V3[mag[,1] == x] }), col = "honeydew", ylab = "difference between counts", main = "per-sample differences", names = trunSampNames, las = 2, varwidth = TRUE, pars = list(cex.axis = 0.8))

boxplot(sapply(as.list(allSampNames),function(x) { mag$V3[mag[,1] == x] }), col = "honeydew", ylab = "difference between counts", main = "per-sample differences", names = trunSampNames, las = 2, ylim = range(-5,5000), varwidth = TRUE, pars = list(cex.axis = 0.8))

boxplot(sapply(as.list(allSampNames),function(x) { mag$V3[mag[,1] == x] }), col = "honeydew", ylab = "difference between counts", main = "per-sample differences", names = trunSampNames, las = 2, ylim = range(-5,20), cex = 0.9, varwidth = TRUE, pars = list(cex.axis = 0.8))
