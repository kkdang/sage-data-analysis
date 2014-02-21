## KKD for Sage Bionetworks
## Sept. 19, 2013
## Looking at QC metrics for common mind project

library("synapseClient")
setwd("~/Computing//commonmind/analysis/")

synapseLogin('kristen.dang@sagebase.org')

#### write out combined file for spotfire
syn_report = synGet('syn2279446')
rna_align_report = read.csv(getFileLocation(syn_report))
syn_file = synGet('syn2279445')
rna_isolation_report = read.csv(getFileLocation(syn_file))

# combine the files
combined_report = cbind(rna_align_report, rna_isolation_report[match(toupper(rna_align_report$Sample.RNA.ID), toupper(rna_isolation_report$Sample.RNA.ID)),4:ncol(rna_isolation_report)])

# test plot
plot(combined_report$RIN, combined_report$Exonic.Rate)

# write file
write.csv(combined_report, file = "QC_metrics_combined.csv")



# Add RIN to other alignment docs
augmented_report = rna_align_report
augmented_report$Mapped.Reads = as.numeric(gsub(",","",augmented_report$Mapped.Reads))
augmented_report$Total.Reads = as.numeric(gsub(",","",augmented_report$Total.Reads))
augmented_report$Exonic.Rate = 0.01*as.numeric(gsub("%","",augmented_report$Exonic.Rate))
augmented_report$rRNA.Rate = 0.01*as.numeric(gsub("%","",augmented_report$rRNA.Rate))
augmented_report$RIN = rna_isolation_report$RIN[match(toupper(augmented_report$Sample.RNA.ID), toupper(rna_isolation_report$Sample.RNA.ID))]

plot(augmented_report$rRNA.Rate, augmented_report$Exonic.Rate*augmented_report$Mapped.Reads/1e6)
plot(augmented_report$rRNA.Rate, augmented_report$Exonic.Rate*augmented_report$Total.Reads/1e6)


# mapped reads vs RIN
plot(augmented_report$RIN, log10(augmented_report$Mapped.Reads), main = "Alignment yield (log10) vs RIN")
plot(augmented_report$RIN, augmented_report$Mapped.Reads, main = "Alignment yield vs RIN")


# fraction mapped vs RIN
augmented_report$Percent.Aligned = as.numeric(substr(augmented_report$Percent.Aligned,1,4))
plot(augmented_report$RIN, augmented_report$Percent.Aligned, main = "Fraction reads aligned vs RIN")


# fraction mapped vs RIN
augmented_report$rRNA.Rate = as.numeric(substr(augmented_report$rRNA.Rate,1,3))
plot(augmented_report$RIN, log10(augmented_report$rRNA.Rate), main = "rRNA aligned vs RIN")

# fraction mapped vs Intragenic
augmented_report$Intragenic.Rate = as.numeric(substr(augmented_report$Intragenic.Rate,1,4))
plot(augmented_report$RIN, augmented_report$Intragenic.Rate, main = "Intragenic aligned vs RIN")

# fraction mapped vs Exonic
augmented_report$Exonic.Rate = as.numeric(substr(augmented_report$Exonic.Rate,1,4))
plot(augmented_report$RIN, augmented_report$Exonic.Rate, main = "Exonic aligned vs RIN")

# fraction mapped vs Intronic
augmented_report$Intronic.Rate = as.numeric(substr(augmented_report$Intronic.Rate,1,4))
plot(augmented_report$RIN, augmented_report$Intronic.Rate, main = "Intronic aligned vs RIN")


# fraction mapped vs Intergenic
augmented_report$Intergenic.Rate = as.numeric(substr(augmented_report$Intergenic.Rate,1,4))
plot(augmented_report$RIN, augmented_report$Intergenic.Rate, main = "Intergenic aligned vs RIN")


# fraction mapped vs Genes detected
augmented_report$Genes.Detected = as.numeric(gsub(",","",augmented_report$Genes.Detected))
plot(augmented_report$RIN, augmented_report$Genes.Detected, main = "Genes detected aligned vs RIN")


# fraction mapped vs Transcripts detected
augmented_report$Transcripts.Detected = as.numeric(gsub(",","",augmented_report$Transcripts.Detected))
plot(augmented_report$RIN, augmented_report$Transcripts.Detected, main = "Transcripts detected aligned vs RIN")


# fraction mapped vs 
augmented_report$Expression.Profiling.Efficiency = as.numeric(substr(augmented_report$Expression.Profiling.Efficiency,1,4))
plot(augmented_report$RIN, augmented_report$Expression.Profiling.Efficiency, main = "Exp profiling efficiency vs RIN")
