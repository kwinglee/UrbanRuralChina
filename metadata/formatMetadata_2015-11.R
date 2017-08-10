##format metadata

rm(list=ls())
setwd("")

md = read.table("metadata_11232015.txt", header=T, sep="\t", colClasses=c(rep("numeric", 49), "character", rep("numeric", 11)))
md = md[,-(1:3)] #remove extra id columns
md = cbind(idind=md$indid, ruralUrban=md$ruralUrban, md[,names(md)!="indid" & names(md) != "ruralUrban"])

##remove missing values
md = md[!is.na(md$hsCRP), ]
##also remove F9, which is all 0s
md = md[,names(md)!="F9"]

##write table
names(md)[1] = "sampleID"
write.table(md, "metadata_cleaned11-2015.txt", sep="\t", row.names=F, col.names=T, quote=F)

##merge with microbiome and metabolome data for classifier
otu = read.table("abundantOTULogNormal.txt", header=T, sep="\t", colClasses=c("character", "character", rep("numeric", 867)))
mb = read.table("MetabolonScaledImpData-transpose.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 337)))
names(mb)[1] = "sampleID"
names(mb)[2] = "ruralUrban"

##metabolite dataset is first timepoint only, so not OTU sample IDs starting with B
##->to get same names, remove leading 0 and ending A from OTU sample names
otu = otu[!grepl("^B", otu$sampleID),]
names = otu$sampleID
names = sub("^0", "", names)
names = sub("A$", "", names)
otu$sampleID = as.numeric(names)

##metadata-microbiome
mer = merge(otu,md, by=c("sampleID","ruralUrban"))
write.table(mer, "OTUWithMetadata_v3.txt", sep="\t", row.names=F, col.names=T, quote=F)

##metadata-metabolome
mer = merge(mb, md, by=c("sampleID","ruralUrban"))
write.table(mer, "MetabolonWithMetadata_v3.txt", sep="\t", row.names=F, col.names=T, quote=F)

##metadata-microbiome-metabolome
mer = merge(otu,md, by=c("sampleID","ruralUrban"))
mer = merge(mer,mb, by=c("sampleID","ruralUrban"))
write.table(mer, "OTUMetadataMetabolon_v3.txt", sep="\t", row.names=F, col.names=T, quote=F)
