# KKD for Sage Bionetworks
# March 27, 2015
# Adds library prep batch to metadata


library('synapseClient')
synapseLogin()

metaEnt = synGet('syn3241418')
meta = read.delim(getFileLocation(metaEnt),header = TRUE)
meta$libraryBatch = rep("NA",times = nrow(meta))

batch1Ent = synGet('syn3387146')
batch1 = read.csv(getFileLocation(batch1Ent),header = TRUE,skip=1)
setdiff(batch1$UDFemale.Investigator.Sample.Name,meta$Px_Code)
intersect(batch1$UDFemale.Investigator.Sample.Name,meta$Px_Code)
meta$libraryBatch[which(meta$Px_Code %in% batch1$UDFemale.Investigator.Sample.Name)] = "batch1"


batch2Ent = synGet('syn3387147')
batch2 = read.csv(getFileLocation(batch2Ent),header = TRUE,skip=1)
setdiff(batch2$UDF.Investigator.Sample.Name,meta$Px_Code)
intersect(batch2$UDF.Investigator.Sample.Name,meta$Px_Code)
meta$libraryBatch[which(meta$Px_Code %in% batch2$UDF.Investigator.Sample.Name)] = "batch2"


batch3Ent = synGet('syn3387148')
batch3 = read.csv(getFileLocation(batch3Ent),header = TRUE,skip=1)
setdiff(batch3$UDF.Investigator.Sample.Name,meta$Px_Code)

x = setdiff(batch3$UDF.Investigator.Sample.Name,meta$Px_Code)[2:6]
y = data.frame(strsplit(x,split = "/"))
match(x,batch3$UDF.Investigator.Sample.Name)

batch3$UDF.Investigator.Sample.Name = as.character(batch3$UDF.Investigator.Sample.Name)
batch3$UDF.Investigator.Sample.Name[match(x,batch3$UDF.Investigator.Sample.Name)] = as.vector(t(y[1,]))

x = grep(pattern = "C1125",x = meta$Px_Code)
y = meta$Px_Code[x]
z = grep(pattern = "C1125",x=batch3$UDF.Investigator.Sample.Name)
batch3$UDF.Investigator.Sample.Name[z] = as.character(y)

batch3$UDF.Investigator.Sample.Name = as.factor(batch3$UDF.Investigator.Sample.Name)
batch3$UDF.Investigator.Sample.Name
setdiff(batch3$UDF.Investigator.Sample.Name,meta$Px_Code)
intersect(batch3$UDF.Investigator.Sample.Name,meta$Px_Code)
meta$libraryBatch[which(meta$Px_Code %in% batch3$UDF.Investigator.Sample.Name)] = "batch3"

write.table(meta,file=getFileLocation(metaEnt),quote = FALSE,row.names = FALSE,sep = "\t")
metaEnt = synStore(metaEnt)
