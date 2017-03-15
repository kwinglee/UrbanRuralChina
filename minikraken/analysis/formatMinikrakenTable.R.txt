##log normalize and add metadata to minikraken results

rm(list=ls())
setwd("")

levels = c("domain", "phylum", "class", "order", "family", "genus", "species")

meta = read.table("genus_taxaAsColumns_relAbunFwd.txt", sep="\t", header=T,
                  colClasses=c(rep("character", 3), rep("numeric", 347)))
meta$sampleID = gsub("_1", "", meta$sampleID)
meta = meta[,1:2]

for(lev in levels) {
  print(lev)
  table = read.table(paste("minikraken_merged_", lev, ".txt", sep=""), 
                     sep="\t", header=F, colClasses="character", quote="") #one species has a ' : Synechococcus_sp._JA-2-3B'a(2-13)
  
  ##remove taxonomy
  table = table[,-2]
  table[1,1] = "sampleID"
  
  ##transpose
  write.table(t(table), paste("minikraken_merged_taxaAsCol_", lev, ".txt", sep=""), 
              sep="\t", row.names = F, col.names = F, quote=F)
  nc = nrow(table)
  table = read.table(paste("minikraken_merged_taxaAsCol_", lev, ".txt", sep=""), 
                     sep = "\t", header=T, colClasses=c("character", rep("numeric", nc-1)), quote="")
  
  ##add number of reads and metadata
  table = merge(meta, table, by="sampleID")
  n = rowSums(table[,3:ncol(table)]) #total number reads in each sample
  table = data.frame(sampleID=table$sampleID, numReads=n, table[,2:ncol(table)], stringsAsFactors = F)
  write.table(table, paste("minikraken_merged_taxaAsCol_", lev, ".txt", sep=""), 
              sep="\t", row.names = F, col.names = T, quote=F)
  
  ##log normalize and write
  lognorm = table
  nc = ncol(lognorm)
  start = 4 #first taxa after metadata
  n = rowSums(table[,start:nc]) #number of reads in each sample
  sumX = sum(n) #total number of reads in all samples = sum(table[,4:nc])
  N = nrow(lognorm) #total number of samples
  for(col in start:nc) {
    for(row in 1:N) {
      lognorm[row, col] = log10(table[row, col]/n[row] * sumX/N + 1)
    }
  }
  
  write.table(lognorm, paste("minikraken_merged_taxaAsCol_logNorm_", lev, ".txt", sep=""),
              sep="\t", row.names = F, col.names = T, quote = F)
}