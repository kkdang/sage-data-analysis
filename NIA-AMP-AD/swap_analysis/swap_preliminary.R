#! /usr/bin/env Rscript
library('synapseClient')
synapseLogin()


library('githubr')
source('/Users/kkdang/Computing/rgithubclient_authenticate.R')
library('edgeR')
library(gplots)
library(RColorBrewer)
library(beeswarm)


sourceRepoFile(sageCode, "rnaseq_analysis_functions.R")

sourceRepoFile(sageCode, "NIA-AMP-AD/swap_analysis/mssm_swap_preliminary.R")
MmayoCountsPalx = mayoCountsPalx
MmssmCountsPalx = mssmCountsPalx
Msamples = c(colnames(MmayoCountsPalx), colnames(MmssmCountsPalx))
rm(mssm, mssmR_counts, mssmPalx.dge, mssmR.dge, mssmCountsPalx)
rm(mayo, mayoR_counts, mayoPalx.dge, mayoR.dge, mayoCountsPalx)
rm(metadata)

sourceRepoFile(sageCode, "NIA-AMP-AD/swap_analysis/rosmap_swap_preliminary.R")
BbroadCountsPalx = broadCountsPalx
BmayoCountsPalx = mayoCountsPalx
BmssmCountsPalx = mssmCountsPalx
Bsamples = c(colnames(BmssmCountsPalx), colnames(BmayoCountsPalx), colnames(BbroadCountsPalx))

rm(broadData, broad, broadCountsPalx, broadPalx.dge, broadR_counts, broadR.dge)
rm(mayoData, mayo, mayoCountsPalx, mayoPalx.dge, mayoR_counts, mayoR.dge)
rm(mssmData, mssm, mssmCountsPalx, mssmPalx.dge, mssmR_counts, mssmR.dge)

sourceRepoFile(sageCode, "NIA-AMP-AD/swap_analysis/mayo_swap_preliminary.R")
YmayoCountsPalx = mayoCountsPalx
YmssmCountsPalx = mssmCountsPalx
Ysamples = c(colnames(YmssmCountsPalx), colnames(YmayoCountsPalx))

rm(mayoData, mayo, mayoCountsPalx, mayoPalx.dge, mayoR_counts, mayoR.dge)
rm(mssmData, mssm, mssmCountsPalx, mssmPalx.dge, mssmR_counts, mssmR.dge)


