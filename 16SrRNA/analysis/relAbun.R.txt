##generate tables of relative abundance

rm(list=ls())

setwd("")

taxaLevels <- c("phylum","class","order","family","genus")

##for given table, returns table converted to relative abundance
relAbun <- function(tab) {
  norm = tab
  nc = ncol(norm)
  start = 4 #first taxa after metadata
  n = rowSums(tab[,start:nc]) #number of reads in each sample
  for(col in start:nc) {
    for(row in 1:nrow(tab)) {
      norm[row, col] = tab[row, col]/n[row]
    }
  }
  return(norm)
}

##for given taxa, returns table of counts with metadata added, forward reads only
getTable <- function(taxa) {
  ##get metadata
  fileName <- paste(taxa, "_taxaAsColumnsLogNorm_WithMetadata.txt", sep ="")
  table <-read.table(fileName,header=TRUE,sep="\t")
  numCols <- ncol(table)
  table <-read.table(fileName,header=TRUE,sep="\t",colClasses=c("character", "numeric", "numeric", "character", "character", rep("numeric", numCols-5)))
  table = table[table$readNumber==1,]
  
  ##get table
  fileName <- paste(taxa, "_taxaAsColumns.txt", sep ="")
  unlog <-read.table(fileName,header=TRUE,sep="\t")
  numCols <- ncol(unlog)
  unlog <-read.table(fileName,header=TRUE,sep="\t",colClasses=c("character", rep("numeric", numCols-1)))
  unlog = unlog[grepl("_1", unlog$sample),]
  names(unlog)[1] = "sampleID"
  meta = data.frame(sampleID = table$sampleID, ruralUrban = table$ruralUrban, timepoint = table$timepoint,
                    stringsAsFactors = F)
  unlog = merge(meta, unlog, by="sampleID")
  
  return(unlog)
}

for(taxa in taxaLevels) {
  print(taxa)
  tab = getTable(taxa)
  norm = relAbun(tab)
  print(rowSums(norm[,-(1:3)]))
  write.table(norm, paste(taxa, "_taxaAsColumns_relAbunFwd.txt", sep=""), sep="\t", row.names=F,
              col.names=T, quote=F)
}

##OTU
table = read.table("abundantOTUForwardTaxaAsColumnsLogNormalWithMetadata.txt", sep="\t", header=T, colClasses=c("character", "numeric", "character", "character", rep("numeric", 870)))
table = cbind(table$sampleID, readNumber=rep(1, nrow(table)), table[,2:ncol(table)]) #add read number
names(table)[1] = "sampleID"
table = table[,-c(6:8)] #remove diversity

unlog = read.table("abundantOTUForwardTaxaAsColumns.txt", sep="\t", header=T, colClasses=c("character", rep("numeric", 867)))
meta = data.frame(sampleID = table$sampleID, ruralUrban = table$ruralUrban, timepoint = table$timepoint,
                  stringsAsFactors = F)
names(unlog)[1] = "sampleID"
unlog = merge(meta, unlog, by="sampleID")

norm = relAbun(unlog)
rowSums(norm[,-(1:3)])
write.table(norm, "abundantOTUForwardTaxaAsColumns_relAbunFwd.txt", sep="\t", row.names=F,
            col.names=T, quote=F)