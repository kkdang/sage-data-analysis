

recounts = read.delim("~/Computing/NIA-AMP-AD/reprocessed/102_120418_filtered_comparison.txt", skip = 1)
head(recounts)
plot(log(recounts$X102_120418.bam), log(recounts$X102_120418_filtered.bam))

re.diff = recounts$X102_120418.bam - recounts$X102_120418_filtered.bam

# How many non-zero differences? = 1647
nonZero = which(abs(re.diff) > 0)
numNonZero =  length(nonZero)

# Out of how many expressed genes? ~0.048 out of 34K expessed
fractionDifferent = numNonZero / length(which(rowSums(recounts[,2:3]) > 0))

# How big are the differences in absolute magnitude?
# Most are < 5, largest is 302
# All diffs are positive, so including chimeric reads (XO:Z:[CUH]T) gives slightly higher counts
max(re.diff[nonZero])
hist(re.diff[nonZero], main = "magnitude of difference in counts", xlab = "difference (in reads)")
hist(re.diff[nonZero], xlim = range(0,50), breaks = 60, main = "magnitude of difference in counts", xlab = "difference (in reads)")


# How big are the differences wrt the expression of the implicated gene?
# Most are less than 0.02 of gene counts
hist(re.diff[nonZero]/ rowMeans(recounts[,2:3])[nonZero] )
hist(re.diff[nonZero]/ rowMeans(recounts[,2:3])[nonZero], xlim = range(0,0.5), breaks = 40 )
hist(re.diff[nonZero]/ rowMeans(recounts[,2:3])[nonZero], xlim = range(0,0.1), breaks = 80 )

# How many have a difference > 10% of gene counts?
length(which(re.diff[nonZero]/ rowMeans(recounts[,2:3])[nonZero] > 0.1))


# How many have a difference > 50% of gene counts?
length(which(re.diff[nonZero]/ rowMeans(recounts[,2:3])[nonZero] > 0.5))

# Where are the ones with large diffs? Some enrichment on chr12 and 18.
namesNonZero = recounts$Geneid[nonZero]
namesNonZero[which(re.diff[nonZero]/ rowMeans(recounts[,2:3])[nonZero] > 0.5)]
write.table(namesNonZero[which(re.diff[nonZero]/ rowMeans(recounts[,2:3])[nonZero] > 0.5)], file = "genes_with_large_count_diffs", quote = FALSE,sep = "\n", row.names = FALSE, col.names = FALSE)
