# KKD for Sage Bionetworks
# Jan. 7, 2014
# QC on GTF - counts file

setwd("~/Computing/commonmind/")

# What genes have novel isoforms? How many have novel isoforms? What is the distribution like? Are ones with a lot of novel isoforms enriched for any pathways?
# First run on command line:
# awk '$3 ~ /j/ {print $1}' 538_merged_with_Ensemble.merged.gtf.tmap | sort | uniq -c > ~/scratch/genes_with_novel_isoforms.v4_merged
genesNovel = read.table("~/Computing/commonmind/data/genes_with_novel_isoforms.v4_merged", row.names = 2)

pdf(file="analysis/hist_novel_isoforms_per_gene.pdf")
hist(as.numeric(genesNovel[,1]), col = "cadetblue2", main = "novel isoforms per gene\nexcluding genes without novel isoforms", xlab = "number")
dev.off()
hist(as.numeric(genesNovel[,1]), col = "cadetblue2", main = "novel isoforms per gene", xlab = "number", xlim = range(0,20), breaks = 160)

library('MASS')
pdf(file="analysis/hist_log_novel_isoforms_per_gene.pdf")
truehist(log10(as.numeric(genesNovel[,1])), col = "cadetblue2", main = "novel isoforms per gene\nexcluding genes without novel isoforms", xlab = "log10(number)")
dev.off()
write.table(rownames(genesNovel)[which(genesNovel[,1] > 50)], file="data/genes_w_gt50_novel.list", quote=FALSE, row.names=FALSE, col.names=FALSE)
# then use this file in CPDB


# Are pseudogenes generating a disproportionate number of novel transcripts?
# First run on command line:
# awk '$3 ~ /j/ {print $2}' 538_merged_with_Ensemble.merged.gtf.tmap | sort | uniq -c > ~/scratch/transcripts_similar_to_novel_isoforms.v4_merged
library(biomaRt)
transNovel = read.table("~/Computing/commonmind/data/transcripts_similar_to_novel_isoforms.v4_merged", row.names = 2)

Hs = useMart("ensembl")
Hs = useDataset("hsapiens_gene_ensembl", Hs)
biotypes = getBM(attributes=c("ensembl_transcript_id","hgnc_symbol","transcript_biotype"),filters="ensembl_transcript_id",values=rownames(transNovel), mart=Hs)

transNovel[,2] = biotypes$transcript_biotype[match(rownames(transNovel),biotypes$ensembl_transcript_id)]
table(transNovel[,2])
pdf(file="analysis/transcripts_similar_to_novel-biotypes.pdf")
dotchart(table(transNovel[,2]),cex = .8, bg = "orangered")
dev.off()