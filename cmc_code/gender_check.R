# KKD for Sage Bionetworks
# Feb. 5, 2014
# Calculates gender based on fraction Y chrom alignments; compares to reported gender

library('synapseClient')
setwd('~/Computing/commonmind/analysis/discord')
synapseLogin()

# get merged metadata file
syn_report = synGet('syn2299154')
merged_report = read.csv(getFileLocation(syn_report))
grep(pattern="Total", x=colnames(merged_report))
colnames(merged_report)[grep(pattern="Total", x=colnames(merged_report))]
head(merged_report)

# Make two input files from the Ychrom_counts.out file: one with counts and one with file names.
# $ awk '$1 !~ /[1-9]/' Ychrom_counts.txt | sed -e 's/samtools view -c \/projects\/CommonMind\/data\/FROM_CORE\/Production\/fastQ\/BySample\/tophat209.2.2\/tophat_//' | sed -e "s/'chrY'//" > Ychrom_names.out
# $ awk '$1 ~ /[1-9]/' Ychrom_counts.txt > Ychrom_data.out

# Read in counts on Y and file names in separate files and combine
Ycounts = read.delim('~/Computing/commonmind/data/Ychrom_data.out', head = F)
hist(log(Ycounts[,1]), main = "distribution of chrY alignments", col = "papayawhip", xlab = "log(Y counts)")
Ynames = read.delim('~/Computing/commonmind/data/Ychrom_names.out', head = F)
rownames(Ycounts) = matrix(unlist(strsplit(as.character(Ynames[,1]), split="/")),ncol = 2, byrow = T)[,1]
Ycounts[,2] = merged_report$DLPFC_RNA_report..Total.Reads[match(rownames(Ycounts), merged_report$DLPFC_RNA_isolation..Sample.RNA.ID)]

# calculated alignment fraction
Ycounts[,3] = log(Ycounts[,1] / Ycounts[,2])

# make sure distribution is bimodal, look for cutoff
hist(Ycounts[,3], col = "papayawhip", xlab = "log(Y counts/total counts)", breaks = 40, main = "distribution of fraction Y counts")
hist(Ycounts[,3], col = "papayawhip", xlab = "log(Y counts/total counts)", breaks = 60, xlim = range(-8,-6))
range(log(Ycounts[,1] / Ycounts[,2]), na.rm = T)

# assign gender
Ycounts[,4] = rep("Unknown", nrow(Ycounts))
males = which(Ycounts[,3] > -7.4)
Ycounts[males,4] = rep("Male", length(males))
females = which(Ycounts[,3] <= -7.4)
Ycounts[females,4] = rep("Female", length(females))
colnames(Ycounts) = c("Y_counts", "total_reads", "Y_fraction", "calc_gender")
Ycounts$calc_gender = as.factor(Ycounts$calc_gender)

# add gender from metadata
grep(pattern="Gender", x=colnames(merged_report))
head(merged_report[,5])
Ycounts[,5] = merged_report$Gender[match(rownames(Ycounts), merged_report$DLPFC_RNA_isolation..Sample.RNA.ID)]
colnames(Ycounts) = c("Y_counts", "total_reads", "Y_fraction", "calc_gender", "reported_gender")
levels(Ycounts$reported_gender) = levels(Ycounts$calc_gender) 

# compare gender assignments
diffGender = Ycounts[,4] != Ycounts[,5]
Ycounts[which(diffGender),]
