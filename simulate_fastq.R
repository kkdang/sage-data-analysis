# Simulating RNAseq reads for workflow comparison
# KKDang for Sage Bionetworks
# July 14, 2015


library(polyester)
library('Biostrings')
setwd('~/Computing/Pipes/')

inTranscripts = '/Users/kristen/Computing/external_data/genomes/gencode.v22.lncRNA-pc_transcripts.fa'
numberOfTranscripts = count_transcripts(inTranscripts)
lengthOfTranscripts = fasta.info(inTranscripts)

# Randomly choose some transcripts to not be expressed
notObserved = sample(x = numberOfTranscripts,size = 0.05*numberOfTranscripts,replace = FALSE)
effectiveLenOfTrans = lengthOfTranscripts
effectiveLenOfTrans[notObserved] = 0

# Getting each transcript's fraction of total bases 
totalReads = 2e5
readLen = 100
transFraction = effectiveLenOfTrans/sum(effectiveLenOfTrans)
hist(log(transFraction))
sum(transFraction)

# Expected reads per transcript
expectedReads = round(transFraction*totalReads)
hist(log(expectedReads))
head(as.numeric(expectedReads))
sum(expectedReads)

length(which(expectedReads == 0))
notSequenced = which(expectedReads == 0)


# subset the FASTA file to only include sequenced transcripts
allFasta = readDNAStringSet(inTranscripts)
sequenced_fasta = allFasta[-notSequenced]
writeXStringSet(sequenced_fasta, paste(getwd(), 'sequencedFasta.fa',sep="/"))
finalExpectedReads = expectedReads[-notSequenced]

simulate_experiment(fasta = 'sequencedFasta.fa',num_reps = c(30,0),reads_per_transcript = finalExpectedReads, fold_changes = rep(1,times=length(finalExpectedReads)),outdir = getwd(),distr='empirical',error_model='illumina5',bias='rnaf',gc_bias=TRUE,size = finalExpectedReads/4)
