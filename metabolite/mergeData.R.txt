##code to generate merged metabolite and OTU tables for classifiers

rm(list = ls())
setwd("")

##get datasets
mb = read.table("MetabolonScaledImpData-transpose.txt", header=T, sep="\t", colClasses=c("character", "character", rep("numeric", 337)))#metabolon
otu = read.table("AbundantOTULogNormal.txt", header=T, sep="\t", colClasses=c("character", "character", rep("numeric", 867)))#otu

###metabolite+OTU
##metabolite dataset is first timepoint only, so not OTU sample IDs starting with B
##->to get same names, remove leading 0 and ending A from OTU sample names
names = otu$sampleID
names = sub("^0", "", names)
names = sub("A$", "", names)
otu$sampleID = names
names(otu)[2] = "classification"
names(mb)[1] = "sampleID"
##check
df1 = data.frame(sampleID=otu$sampleID, ruralUrban=otu$classification)
df2 = data.frame(sampleID=mb$sampleID, class=mb$classification)
mer = merge(df1, df2, by="sampleID")
sum(mer$ruralUrban!=mer$class)
sum(mer$ruralUrban==mer$class)
##merge columns
mer = merge(mb, otu, by=c("sampleID","classification"))
write.table(mer, "MetabolonWithOTU.txt", sep="\t", row.names=F, col.names=T, quote=F)
