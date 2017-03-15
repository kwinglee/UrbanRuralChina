##supplemental table of number of taxa identified at each level for each technique

rm(list=ls())
setwd("")

results = data.frame(level = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU/Species", "KEGG Module", "KEGG Pathway"),
                     rrna = rep(NA, 9),
                     wgs = rep(NA, 9))

####16S
taxaLevels = c("phylum","class","order","family","genus", "otu")
for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  print(taxa)
  inFileName = paste(taxa,  "_taxaAsColumnsLogNorm_WithMetadata.txt", sep ="")
  if(taxa=="otu") {
    inFileName = "abundantOTUForwardTaxaAsColumnsLogNormalWithMetadata.txt"
  }
  table = read.table(inFileName,header=TRUE,sep="\t")
  numCols = ncol(table)
  myColClasses = c(rep("character",5), rep("numeric", numCols-5))
  if(taxa=="otu") { #missing read number
    myColClasses = c(rep("character",4), rep("numeric", numCols-4))
  }
  table = read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)
  
  ##add read number to OTU and remove diversity measures
  if(taxa=="otu") {
    ##add read number
    readNumber = rep("1", nrow(table))
    table = cbind(table[,1], readNumber, table[,-1])
    names(table)[1] = "sampleID"
    ##remove diversity measures
    table = table[,-(6:8)]
    numCols = ncol(table)
  }
  
  ##forward reads only
  table = table[table$readNumber=="1",] #some taxa pass next test when include all reads
  
  ##filter columns > 1/4 samples are 0
  filter <- function(c) (return(sum(table[,c] != 0) > nrow(table) / 4))
  keep = sapply(1:ncol(table), filter)
  table = table[,keep]
  
  results$rrna[i+1] = ncol(table)-5
}

####minikraken
taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")
for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  print(taxa)
  file = paste("minikraken_merged_taxaAsCol_logNorm_", taxa, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  ##filter zeros in > 1/4 samples
  filter <- function(i) (return(sum(table[,i] != 0 ) > nrow(table) / 4))
  keep = sapply(1:ncol(table), filter)
  table = table[,keep]
  
  results$wgs[i] = ncol(table)-3
}

####humann
levels = c("module", "pathway")
for(i in 1:length(levels)) {
  lev = levels[i]
  print(lev)
  file = paste("humann_keggAsCol_withRurUrb_", lev, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", ncol-2)))
  
  ##filter zeros in > 1/4 samples
  filter <- function(i) (return(sum(table[,i] != 0 ) > nrow(table) / 4))
  keep = sapply(1:ncol(table), filter)
  table = table[,keep]
  
  results$wgs[i+length(taxaLevels)] = ncol(table)-6
}

names(results) = c("level", "16S Sequencing", "Whole Genome Sequencing")
write.table(results, "numTaxa.txt", sep="\t", row.names = F, col.names = T, quote = F)