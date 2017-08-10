##generate supplemental table of sequence stats and accession numbers
##columns: sample name, subject ID, urban/rural, timepoint, SRA 16S accession, SRA WGS accession,
##    forward MG-RAST 16S ID, reverse 16S ID, number of reads 16S, number reads WGS
##last rows are average and standard deviation

rm(list=ls())
setwd("")

##function to clean sample IDs
cleanSampleID <- function(sID) {
  ###clean up sample ids
  ##for time 2, remove A at end, move B there
  sID[grepl("B", sID)] = sub("A", "", sID[grepl("B", sID)])
  sID[grepl("B", sID)] = paste(sub("B", "", sID[grepl("B", sID)]), "B", sep="")
  
  ##add leading 0s
  sID[nchar(sID) == 2] = paste("0", sID[nchar(sID) == 2], sep="")
  sID[nchar(sID) == 3] = paste("0", sID[nchar(sID) == 3], sep="")
  return(sID)
}

#####get metadata
meta = read.table("phylum_taxaAsColumnsLogNorm_WithMetadata.txt",
                  header=T, sep="\t", stringsAsFactors = F)
meta = meta[grepl("_1", meta$sampleID),]
meta = data.frame(sampleID = cleanSampleID(sub("_1", "", meta$sampleID)),
                  subjectID = meta$patientID,
                  ruralUrban = meta$ruralUrban,
                  timepoint = ifelse(meta$timepoint=="first_A", "first", "second"),
                  stringsAsFactors = F)


#####SRA accession numbers
sra16s = read.table("China_16S_srametadata.txt-processed-ok.tsv", header=T, sep="\t")
sraWGS = read.table("China_wgs_srametadata.txt-processed-ok.tsv", header=T, sep="\t")

sra16s = data.frame(sampleID = sra16s$sample_name, 
                    accession16s = as.character(sra16s$accession),
                    stringsAsFactors = F)
sraWGS = data.frame(sampleID = sub("_wgs", "", sraWGS$sample_name), 
                    accessionWGS = as.character(sraWGS$accession),
                    stringsAsFactors = F)

##add second timepoint to WGS with NA for accession number
sraWGS = rbind(sraWGS, 
               data.frame(sampleID = sra16s$sampleID[!(sra16s$sampleID %in% sraWGS$sampleID)],
                          accessionWGS = NA,
                          stringsAsFactors = F))

#####MGRAST columns
####get MG-RAST ID
mgTab = read.table("2016-2-10table_all160_from_mgrast.txt", sep="\t", header=T, 
                   colClasses=c("character", "character", "numeric", "numeric", rep("character", 6)))
mgTab = data.frame(file=mgTab$Metagenome.Name, 
                   MG.RAST.ID=mgTab$MG.RAST.ID, 
                   # raw.counts=mgTab$Sequence.Count, 
                   stringsAsFactors = F)

####convert MG-RAST to sample IDs; make forward and reverse into separate columns; remove file name
mg1 = mgTab[grepl("_forward", mgTab$file),]
mg2 = mgTab[grepl("_reverse", mgTab$file),]
names(mg1) = c("sampleID", "MG.RAST.FWD")
names(mg2) = c("sampleID", "MG.RAST.REV")
mg1$sampleID = sub("_forward", "", mg1$sampleID)
mg2$sampleID = sub("_reverse", "", mg2$sampleID)

####get number of 16S reads
readTable = read.table("number16Sreads.txt",
                       header=T, sep="\t", stringsAsFactors = F)
readTable$file = sub("first_", "", sub("second_", "", sub(".fq.gz", "", readTable$file)))
##check forward and reverse have same number reads
r1 = readTable[grepl("_1", readTable$file),]
r2 = readTable[grepl("_2", readTable$file),]
r1$file = sub("_1", "", r1$file)
r2$file = sub("_2", "", r2$file)
readMrg = merge(r1, r2, by="file")
any(readMrg$numberReads.x != readMrg$numberReads.y) #FALSE -> pass
##format sample ids
readTable = data.frame(sampleID = cleanSampleID(r1$file),
                       number16sReads = r1$numberReads,
                       stringsAsFactors = F)

#####add number of WGS reads
readTable$numberWGSreads = ifelse(grepl("A", readTable$sampleID), 10000000, NA)

#####merge
all = merge(meta, sra16s, by="sampleID")
all = merge(all, sraWGS, by="sampleID")
all = merge(all, mg1, by="sampleID")
all = merge(all, mg2, by="sampleID")
all = merge(all, readTable, by="sampleID")

all = all[order(all$sampleID),]

#####add average and standard deviation
avg = c("Average", rep("", ncol(all)-1))
sd = c("Standard deviation", rep("", ncol(all)-1))
for(i in 9:ncol(all)) {
  avg[i] = mean(all[,i], na.rm=T)
  sd[i] = sd(all[,i], na.rm=T)
}
all = rbind(all, avg, sd)

write.table(all, "seqStats.txt", row.names=F, col.names=T, quote=F, sep="\t")
