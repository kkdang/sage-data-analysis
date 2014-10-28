#! /bin/bash
# Making VCF tracks for cranio dataset on belltown

# install current bcftools
wget http://sourceforge.net/projects/samtools/files/samtools/1.1/bcftools-1.1.tar.bz2/download
mv download bcftools-1.1.tar.bz2
tar -xvjf bcftools-1.1.tar.bz2 
cd bcftools-1.1/

## Make lists of files for cases and controls
ls *vcf.gz | grep -f ../../../../Data/sampleList_cases.txt > ../../../../Data/fileList_cases_vcf.txt 
ls *vcf.gz | grep -f ../../../../Data/sampleList_controls.txt > ../../../../Data/fileList_controls_vcf.txt 


## make separate VCFs for cases and controls
nohup ~/Software/bcftools-1.1/bcftools merge -m both -o ../../../../Data/controls_newmerge.vcf.gz -O z -l ../../../../Data/fileList_controls_vcf.txt &> ../../../../Data/newmerge_contols.log &

nohup ~/Software/bcftools-1.1/bcftools merge -m both -o ../../../../Data/cases_newmerge.vcf.gz -O z -l ../../../../Data/fileList_cases_vcf.txt &> ../../../../Data/newmerge_cases.log &
 
## index
~/Software/bcftools-1.1/bcftools index -t cases_newmerge.vcf.gz 
~/Software/bcftools-1.1/bcftools index -t controls_newmerge.vcf.gz 

## Get any variants found in cases not controls
nohup ~/Software/bcftools-1.1/bcftools isec -O z -p test_comp/ -c both -C cases_newmerge.vcf.gz controls_newmerge.vcf.gz 


## Try using bedtools 
nohup cat controls_merged.vcf > ~/Software/vcftools_0.1.11/bin/vcf-sort > controls_merged_sorted.vcf &
nohup cat cases_merged.vcf > ~/Software/vcftools_0.1.11/bin/vcf-sort > cases_merged_sorted.vcf &
nohup bedtools intersect -v -sorted -a cases_merged_sorted.vcf -b controls_merged_sorted.vcf > ../bedtools_isec/test_sorted.out &


## Data from dbSNP
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/All.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/All.vcf.gz.tbi


## Comparing to dbSNP
nohup ~/Software/bcftools-1.1/bcftools isec -O z -p case-not-dbsnp/ -c none -C cases_merged.vcf.gz All.vcf.gz &
nohup ~/Software/bcftools-1.1/bcftools isec -O z -p case-not-dbsnp-both/ -c both -C cases_merged.vcf.gz All.vcf.gz &