##format RDP results for classifiers

rm(list=ls())
setwd("")

group = c("phylum", "class", "order", "family", "genus")
for(g in group) {
  table = read.table(paste(g, "_taxaAsColumnsLogNorm_WithMetadata.txt", sep=""), sep="\t", header=T)
  ##only keep forward reads
  samps = table$sampleID
  rows = rep(T, length(samps))
  rows[grepl("_2$", samps, fixed=F)] = F
  out = table[rows,]
  out = out[,-c(2,3,5)] #remove all metadata but ruralUrban
  write.table(out, paste(g, "_forSVMLight.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}