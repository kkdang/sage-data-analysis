#! /bin/bash
# Oct. 8, 2014
# Kristen K. Dang for Sage Bionetworks
# Make bigWig tracks via bedGraph files

# Convert GTF to bed, separately for forward and reverse strands
cut -f1 -d";" /sc/orga/projects/GCF/REFERENCES/hg19/ensemble/Homo_sapiens.GRCh37.70.processed.gtf | sed -e 's/gene_id//g' | awk '$7 ~ /-/' |  sort -k1,1 -k4,4n | bedtools merge -s -c 9 -o distinct -i - > Homo_sapiens.GRCh37.70.processed_merged_minusstrand_distinct.bed

cut -f1 -d";" /sc/orga/projects/GCF/REFERENCES/hg19/ensemble/Homo_sapiens.GRCh37.70.processed.gtf | sed -e 's/gene_id//g' | awk '$7 ~ /+/' |  sort -k1,1 -k4,4n | bedtools merge -s -c 9 -o distinct -i - > Homo_sapiens.GRCh37.70.processed_merged_plusstrand_distinct.bed


# Altnerative -- stranded merge of both strands in one file
#cut -f1 -d";" Homo_sapiens.GRCh37.70.processed.gtf | sed -e 's/gene_id//g' | sort -k1,1 -k4,4n | bedtools merge -s -c 9 -o distinct -i - > Homo_sapiens.GRCh37.70.processed_merged_stranded_distinct.bed


module load python py_packages
synapse get syn2713767

# Get list of genes in input (data) file
cut -f1 DLPFC.ensembl.DxSCZ.DE.KNOWN_AND_SVA.ADJUSTED.VOOM_NORMALIZED.LIMMA_FIT.tsv | sort | uniq > genes_in_normfile

# Are there duplicate genes in input (data) file?
cut -f1 DLPFC.ensembl.DxSCZ.DE.KNOWN_AND_SVA.ADJUSTED.VOOM_NORMALIZED.LIMMA_FIT.tsv | sort | uniq -d | head

# Reduce the bed file to contain only the genes in the input file
bsub -q scavenger -W 02:00 -o longgrep.stdout -e longgrep.stderr -L /bin/bash ./long_grep.sh


# Make chrom sizes file required to convert to wig
# *** Maybe change to remove "chr" from chromosome names?
samtools view -H /sc/orga/projects/CommonMind/data/FROM_CORE/Production/fastQ/BySample/tophat209.2.2/tophat_MSSM_RNA_PFC_347/MSSM_RNA_PFC_347.sort.bam > header.txt
cut -f2,3 header.txt | sed -e 's/[SN:|LN:]//g' > chromSizes.txt




# Add the DE data to the bed file, creating separate tracks.
awk 'NF == 4' Homo_sapiens.GRCh37.70.processed_minusstrand_reduced.bed | sed -e 's/[",]//g' > Homo_sapiens.GRCh37.70.processed_minusstrand_reduced_nodual.bed
awk 'NF == 4' Homo_sapiens.GRCh37.70.processed_plusstrand_reduced.bed | sed -e 's/[",]//g' > Homo_sapiens.GRCh37.70.processed_plusstrand_reduced_nodual.bed

# run in own directory, separate from below job
bsub -q scavenger -W 00:10 ~/bin/add_values_to_bedGraph.py DLPFC.ensembl.DxSCZ.DE.KNOWN_AND_SVA.ADJUSTED.VOOM_NORMALIZED.LIMMA_FIT.tsv Homo_sapiens.GRCh37.70.processed_minusstrand_reduced_nodual.bed

bsub -q scavenger -W 00:10 ~/bin/add_values_to_bedGraph.py DLPFC.ensembl.DxSCZ.DE.KNOWN_AND_SVA.ADJUSTED.VOOM_NORMALIZED.LIMMA_FIT.tsv Homo_sapiens.GRCh37.70.processed_plusstrand_reduced_nodual.bed



# Convert to bigwig
module load ucsc-utils
cut -f1,2,3,4 logFCTrack.bed | sort -k1,1 -k2,2n > logFCTrack_sorted.bed
bedGraphToBigWig logFCTrack_sorted.bed chromSizes.txt logFCTrack.bigwig

cut -f1,2,3,4 minusstrand/aveExpTrack.bed | sort -k1,1 -k2,2n > minusstrand/aveExpTrack_sorted.bed
bedGraphToBigWig minusstrand/aveExpTrack_sorted.bed chromSizes.txt minusstrand/aveExpTrack.bigwig

cut -f1,2,3,4 minusstrand/pvalTrack.bed | sort -k1,1 -k2,2n > minusstrand/pvalTrack_sorted.bed         
bedGraphToBigWig minusstrand/pvalTrack_sorted.bed chromSizes.txt minusstrand/pvalTrack.bigwig


cut -f1,2,3,4 plusstrand/logFCTrack.bed | sort -k1,1 -k2,2n > plusstrand/logFCTrack_sorted.bed
bedGraphToBigWig plusstrand/logFCTrack_sorted.bed chromSizes.txt plusstrand/logFCTrack.bigwig

cut -f1,2,3,4 plusstrand/aveExpTrack.bed | sort -k1,1 -k2,2n > plusstrand/aveExpTrack_sorted.bed
bedGraphToBigWig plusstrand/aveExpTrack_sorted.bed chromSizes.txt plusstrand/aveExpTrack.bigwig

cut -f1,2,3,4 plusstrand/pvalTrack.bed | sort -k1,1 -k2,2n > plusstrand/pvalTrack_sorted.bed         
bedGraphToBigWig plusstrand/pvalTrack_sorted.bed chromSizes.txt plusstrand/pvalTrack.bigwig

