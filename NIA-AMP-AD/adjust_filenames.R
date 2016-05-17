mssm = synGet('syn5553119')
mssmData = read.table(getFileLocation(mssm), row.names = 1)
colnames(mssmData)


key = read.csv("~/Downloads/111615_Key_SampleSwap_MayoSamples_ToMtSinai.csv")
tmp = gsub("X","",colnames(mssmData))
tmp2 = gsub("_", "-", tmp)
adjColnames = gsub("[.]", "-", tmp2)
newColnames = key$New.ID[match(adjColnames,as.character(key$Old.ID))]
colnames(mssmData) = newColnames
