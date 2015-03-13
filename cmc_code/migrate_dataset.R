# Removing samples from full dataset


library('synapseClient')
synapseLogin()

metadataEnt = synGet('syn3346482')
metadata = read.csv(getFileLocation(metadataEnt))
colnames(metadata)

filesToMigrate = c('syn2340130', 'syn2340132','syn2342352','syn2338623')

act = Activity(name="Align and Quantitate Reads",
                used=list(
                  list(entity='syn3346718', wasExecuted=T),
                  list(entity='syn3280463', wasExecuted=F)))
act = storeEntity(act)


for (i in 1:length(filesToMigrate)) {
  dataEntity = synGet(filesToMigrate[i])
  dataFile  = read.delim(dataEntity@filePath, header=T, row.names = 1)
  
  ncol(dataFile)
  dataFile_release = dataFile[,which(colnames(dataFile) %in% metadata$DLPFC_RNA_Sequencing_Sample_ID)]
  ncol(dataFile_release)
  head(colnames(dataFile_release))
  tail(colnames(dataFile_release))
  
  nameOptions = c("geneExpressionRaw", "geneExpressionRawUCSC","exonExpressionRaw","geneExpressionNorm")
  fileBase = 'CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_'
  outFileName = paste(fileBase, nameOptions[i],".tsv",sep = "")
  write.table(x = dataFile_release,file = outFileName, quote = FALSE,row.names = TRUE, col.names = TRUE,sep = "\t")
  library(R.utils)
  gzip(outFileName)
  fileEntity = File(path=paste(outFileName, "gz", sep = "."), parentId='syn3275215')

  
  annotList = list(
    geneList=list(normalizationStatus=FALSE, summaryLevel="gene", geneIdentifiers="Ensembl",consortium='CMC', Cohort='MSSM-Penn-Pitt', dataType='mRNA', disease=c('Control','Schizophrenia', 'Bipolar Disorder'), platform='IlluminaHiSeq2500', dissectionSource='CMC DLPFC Dissection 1', SampleSize='613',tissueType='Dorsolateral Prefrontal Cortex', tissueAbreviation='DLPFC', fileType='tsv', summaryType='counts', organism='Homo sapiens', dataContact='info@commonmind.org', migratedSynIDs=filesToMigrate[i], migratedSynVersion=paste('Version',dataEntity$properties$versionNumber,sep = ' '), migratedScript='url', Release='1.0'),
    ucscList=list(normalizationStatus=FALSE, summaryLevel="gene", geneIdentifiers="HGNC",consortium='CMC', Cohort='MSSM-Penn-Pitt', dataType='mRNA', disease=c('Control','Schizophrenia', 'Bipolar Disorder'), platform='IlluminaHiSeq2500', dissectionSource='CMC DLPFC Dissection 1', SampleSize='613',tissueType='Dorsolateral Prefrontal Cortex', tissueAbreviation='DLPFC', fileType='tsv', summaryType='counts', organism='Homo sapiens', dataContact='info@commonmind.org', migratedSynIDs=filesToMigrate[i], migratedSynVersion=paste('Version',dataEntity$properties$versionNumber,sep = ' '), migratedScript='url', Release='1.0'),
    exonList=list(normalizationStatus=FALSE, summaryLevel="exon", geneIdentifiers="Ensembl",consortium='CMC', Cohort='MSSM-Penn-Pitt', dataType='mRNA', disease=c('Control','Schizophrenia', 'Bipolar Disorder'), platform='IlluminaHiSeq2500', dissectionSource='CMC DLPFC Dissection 1', SampleSize='613',tissueType='Dorsolateral Prefrontal Cortex', tissueAbreviation='DLPFC', fileType='tsv', summaryType='counts', organism='Homo sapiens', dataContact='info@commonmind.org', migratedSynIDs=filesToMigrate[i], migratedSynVersion=paste('Version',dataEntity$properties$versionNumber,sep = ' '), migratedScript='url', Release='1.0'),
    fpkmList=list(normalizationStatus=TRUE, summaryLevel="gene", geneIdentifiers="Ensembl", normalizationType='FPKM',consortium='CMC', Cohort='MSSM-Penn-Pitt', dataType='mRNA', disease=c('Control','Schizophrenia', 'Bipolar Disorder'), platform='IlluminaHiSeq2500', dissectionSource='CMC DLPFC Dissection 1', SampleSize='613',tissueType='Dorsolateral Prefrontal Cortex', tissueAbreviation='DLPFC', fileType='tsv', summaryType='counts', organism='Homo sapiens', dataContact='info@commonmind.org', migratedSynIDs=filesToMigrate[i], migratedSynVersion=paste('Version',dataEntity$properties$versionNumber,sep = ' '), migratedScript='url', Release='1.0')
    )
                   
  synSetAnnotations(fileEntity) = annotList[[i]]
  generatedBy(fileEntity) = act
  fileEntity = synStore(fileEntity)
  
}

#x = list(consortium='CMC', Cohort='MSSM-Penn-Pitt', dataType='mRNA', disease=c('Control','Schizophrenia', 'Bipolar Disorder'), platform='IlluminaHiSeq2500', dissectionSource='CMC DLPFC Dissection 1', SampleSize='613',tissueType='Dorsolateral Prefrontal Cortex', tissueAbreviation='DLPFC', fileType='tsv', summaryType='counts', organism='Homo sapiens', dataContact='info@commonmind.org', migratedSynIDs=filesToMigrate[i], migratedSynVersion=paste('Version',dataEntity$properties$versionNumber,sep = ' '), migratedScript='url', Release='1.0')
