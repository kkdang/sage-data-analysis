# Feb. 29, 2016
# KKD for Sage Bionetworks
# Calculates prevalence of each input mutation in sample groups from unsupervised clustering results

setwd('~/Computing/cranio/variants/')

library(synapseClient)
synapseLogin()

library('rGithubClient')
source('/Users/kristen/Computing/external_software/rgithubclient_authenticate.R')

library(biomaRt)


# get sample groups from unsupervised gene-level clustering
sgEnt = synGet('syn4909512')
geneSampleGroups = read.delim(getFileLocation(sgEnt),header = FALSE)
head(geneSampleGroups)
table(geneSampleGroups[,2])


# get key for seq id <-> sample id
sourceRepoFile(sageCode, "scri_cran/process_metadata.R")
metadataFiltered = processMetadata(ver = 5,plot = FALSE)
head(metadataFiltered)
geneSampleGroups[,3] = metadataFiltered$SAMPLE[match(geneSampleGroups[,1],metadataFiltered$Px_Code)]
geneSampleGroups[,2] = as.factor(geneSampleGroups[,2])

sampleGroups = list()
sampleGroups$`2` = geneSampleGroups[,3][geneSampleGroups[,2] == "2"]
sampleGroups$`3` = geneSampleGroups[,3][geneSampleGroups[,2] == "3"]
sampleGroups$`4` = geneSampleGroups[,3][geneSampleGroups[,2] == "4"]
sampleGroups$`5` = geneSampleGroups[,3][geneSampleGroups[,2] == "5"]
sampleGroups$`9` = geneSampleGroups[,3][geneSampleGroups[,2] == "9"]
sampleGroups$`10` = geneSampleGroups[,3][geneSampleGroups[,2] == "10"]
sampleGroups$`13` = geneSampleGroups[,3][geneSampleGroups[,2] == "13"]
sampleGroups$`14` = geneSampleGroups[,3][geneSampleGroups[,2] == "14"]
sampleGroups$`other` = as.vector(na.omit(setdiff(metadataFiltered$SAMPLE, geneSampleGroups[,3])))


sampleNameOrderFromVCF = c("95445", "95446", "95447", "95448", "95450", "95452", "95453", "95454", "95456", "95457", "95458", "95459", "95460", "95462", "95463", "95464", "95466", "95467", "95468", "95469", "95470", "95471", "95472", "95474", "95475", "95476", "95477", "95479", "95481", "95482", "95484", "95485", "95486", "95487", "95488", "95490", "95491", "95492", "95493", "95494", "95495", "95496", "95498", "95500", "95501", "95502", "95505", "95506", "95507", "95509", "95510", "95511", "95512", "95514", "95515", "95516", "95518", "95519", "95520", "95521", "95522", "95524", "95525", "95526", "95527", "95528", "95529", "95530", "95531", "95534", "95535", "95536", "95537", "95540", "95542", "95543", "95544", "95545", "95546", "95547", "95549", "95550", "95551", "95552", "95555", "95556", "95557", "95558", "95559", "95560", "95562", "95563", "95564", "95565", "95566", "95567", "95569", "95570", "95571", "95572", "95573", "95575", "95576", "95578", "95579", "95580", "95581", "95583", "95584", "95585", "95586", "95587", "95588", "95589", "95591", "95592", "95593", "95594", "95595", "95596", "95597", "95598", "95599", "95600", "95601", "95602", "95603", "95604", "95605", "95606", "95609", "95610", "95611", "95612", "95614", "95615", "95616", "95617", "95619", "95620", "95622", "95623", "95624", "95626", "95627", "95628", "95629", "95630", "95631", "95633", "95634", "95635", "95636", "95637", "95638", "95639", "95640", "95642", "95644", "95645", "95646", "95649", "95650", "95651", "95652", "95653", "95654", "95655", "95656", "95657", "95658", "95661", "95662", "95663", "95665", "95668", "95669", "95670", "95671", "95673", "95675", "95676", "95677", "95678", "95679", "95680", "95681", "95682", "95683", "95684", "95685", "95686", "95688", "95689", "95690", "95691", "95693", "95694", "95696", "95697", "95698", "95699", "95700", "95702", "95703")



# Get variant prevalence
exacEnt = synGet('syn5698384')
exac = read.delim(getFileLocation(exacEnt),header = FALSE)
head(exac)
colnames(exac)
colnames(exac)[4:208] = sampleNameOrderFromVCF

groupPrev = apply(exac[,4:208],MARGIN = 1,function(x){calcGroupPrevalence(gt_vec = x)})
head(groupPrev)
hist(groupPrev[1,])
ofInterest = which(groupPrev[1,] < 0.01)
sigDistributions = cbind(exac[ofInterest,1:3], geneNames[ofInterest,2], t(groupPrev[2:10,ofInterest]))
colnames(sigDistributions)[1:4] = c("chrom", "pos", "type", "symbol")

# Get gene symbols for relevant SNPs
Hs=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
Hs = useDataset("hsapiens_gene_ensembl", Hs)
geneNames = matrix(NA, ncol = 2, nrow = nrow(exac))
for (i in 1:nrow(exac)){ 
  x = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),filters=c("chromosome_name", "start", "end"),values=as.list(c(as.character(exac[i,1]), exac[i,2], exac[i,2])), mart=Hs)
  if (nrow(x) == 0) { next }
  if (nrow(x) > 1) {
    geneNames[i,1] = paste(x$ensembl_gene_id, collapse = "-")
    geneNames[i,2] = paste(x$hgnc_symbol, collapse = "-")
  }
  else {
    geneNames[i,1] = x$ensembl_gene_id
    geneNames[i,2] = x$hgnc_symbol
  }  
}

write.table(x = sigDistributions, file = "variants_vs_cluster_groups.txt",append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
