#! /usr/bin/env Rscript
# Dec. 29, 2015
# KKD for Sage Bionetworks


library(synapseClient)
library('githubr')
source('/home/ubuntu/rgithubclient_authenticate.R')
sageCode = getRepo(repository="kkdang/sage-data-analysis")
library(metatest)
library(plyr)
library(doParallel)
setwd('/mnt/')
synapseLogin()
library(argparse)

parser = ArgumentParser(description='Collect and sum assay level counts into a matrix file.')
parser$add_argument('--input', type="character", required=TRUE, help='Input data file.')
parser$add_argument('--parallelJobs', default = 1, type="integer", help='Number of jobs to run concurrently [default %(default)d].')
parser$add_argument('--verbose', default = FALSE, type="binary", help='Whether to print detailed output, at expense of overall runtime. [default %(default)s].')
parser$add_argument('--wd', default = getwd(), type="character", help='Directory for output files [default %(default)s].')
args = parser$parse_args()


## Filter and convert metadata
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)
metadataFiltered$Sample_Type = as.character(metadataFiltered$Sample_Type)
metadataFiltered$Sample_Type[grep("Coronal",metadataFiltered$Sample_Type)] = "Coronal"
metadataFiltered$Sample_Type[grep("Metopic",metadataFiltered$Sample_Type)] = "Metopic"
metadataFiltered$Sample_Type[grep("Sagittal",metadataFiltered$Sample_Type)] = "Sagittal"
metadataFiltered$Sample_Type = as.factor(metadataFiltered$Sample_Type)

## Remove outliers
outlierTable = synTableQuery('SELECT * FROM syn3354503')
outlierData = outlierTable@values
noOutliers.meta = metadataFiltered[-which(metadataFiltered$Px_Code %in% outlierData$transcriptsSet),]


## Restrict to model variables, set column names and order as required by meta-diff
modelToKeep = c("Age_mos.","PCT_CORRECT_STRAND_READS","Initial_growth_duration_days", "SAMPLE", "caseStatus", "Px_Code")
colnames(noOutliers.meta)
model.meta = noOutliers.meta[,which(colnames(noOutliers.meta) %in% modelToKeep)]
head(model.meta)
model.meta$"C_group" = rep(0,nrow(model.meta))
model.meta$C_group[which(model.meta$caseStatus == "case")] = 1
tail(model.meta)

model.meta$Age_mos. = formatC(output.file$Age_mos., digits=6, format="fg")
model.meta$Initial_growth_duration_days = formatC(output.file$Initial_growth_duration_days, digits=6, format="fg")
model.meta$Age_mos. = formatC(output.file$Age_mos., digits=6, format="fg")
head(model.meta)


table(metadataFiltered$Sample_Type)
coronalSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Coronal"])
controlSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Control"])
lambdoidSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Lambdoid"])
metopicSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Metopic"])
sagittalSampleIds = as.character(metadataFiltered$Px_Code[metadataFiltered$Sample_Type == "Sagittal"])

####

fpkm = read.delim("/mnt/allFPKM.dat", header = FALSE, sep = "\t")
colnames(fpkm) = c("transcript_id", "y", "conf-lo","conf-hi","status","SAMPLE")

nodes = 4

modelF = as.formula('y~factor(C_GROUP)+AGE_MOS.+PCT_CORRECT_STRAND_READS+INITIAL_GROWTH_DURATION_DAYS')


run_metatest<-function(mat,metaData,model=modelF){
  tryCatch({  
	z = join(mat,metaData,by="SAMPLE",type = "inner")
	if (sum(z$y == 0) { return(NULL) }
    meta_obj <- metatest(formula = model, variance = variance, data=z)
    res_vec <- c(meta_obj$convergence, meta_obj$tval, meta_obj$dfttest, meta_obj$pttest, meta_obj$bartLLR, meta_obj$pBartlett)
    cov_name <- dimnames(meta_obj$coefficients)[[1]]
    names_vec <- c('Convergence', paste('tval_',cov_name,sep=''), 'dfttest', paste('pttest_',cov_name,sep=''), paste('bartLLR_',cov_name,sep=''), paste('pBartlett_',cov_name,sep=''))
    names(res_vec)<-names_vec
    return(res_vec)}, error = function(e) {
      warning(paste('Metatest failed at Feature: ', mat[1,1],sep=''))
      warning(e)
      return(NULL)
    }
  )
}

doParallel::registerDoParallel(cores = args$parallelJobs)
df_res <- ddply(fpkm, .(Feature), .fun = run_metatest, .parallel = T, .inform = args$inform, .paropts = list(.packages='metatest'))

write.table(df_res, file = "metatest_preprocessed.txt", row.names = F, quote=F, sep='\t')




