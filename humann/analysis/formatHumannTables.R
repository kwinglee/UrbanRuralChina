##format and log humann results
##add metadata

rm(list=ls())

setwd("")

humNames = c("mpm", "mpt")
names = c("module", "pathway")

for(i in 1:length(humNames)) {
  print(humNames[i])
  table = read.table(paste("04b-hit-keg-", humNames[i], "-cop-nul-nve-nve.txt", sep=""), sep="\t", 
                     stringsAsFactors = F, quote = "") #without quote mpt throws errors at ko05016 (Huntington's disease)
  
  ##fix sample names
  table[1,1] = "sampleID" #ID causes errors in Excel
  table[1,] = gsub("kegg_merge_filter_human_", "", gsub(".hit.keg.mp[t|m].cop.nul.nve.nve", "", table[1,]))
  
  ##remove empty rows
  if(table$V1[15] != "InverseSimpson") {
    print(paste("check rows:", h))
  }
  table = table[-(2:14),]
  
  ##write so can open in Excel, and transpose
  write.table(table, paste("humann_keggAsRow_", names[i], ".txt", sep=""), sep="\t", row.names = F, 
              col.names = F, quote=F)
  trans = t(table)
  write.table(t(table), paste("humann_keggAsCol_allWithName_", names[i], ".txt", sep=""), sep="\t", row.names = F, 
              col.names = F, quote=F)
  write.table(trans[-2,], paste("humann_keggAsCol_", names[i], ".txt", sep=""), sep="\t", row.names = F, 
              col.names = F, quote=F) #second row with descriptions will get in the way later
  
  ##add metadata
  nr = nrow(table)
  table = read.table(paste("humann_keggAsCol_", names[i], ".txt", sep=""), sep="\t", header=T,
                     colClasses = c("character", rep("numeric", nr-1)))
  meta = read.table("genus_taxaAsColumns_relAbunFwd.txt", sep="\t", header=T,
                    colClasses=c(rep("character", 3), rep("numeric", 347)))
  meta = meta[meta$timepoint=="first_A",1:2]
  meta$sampleID = gsub("_1", "", meta$sampleID)
  table = merge(meta, table, by="sampleID")
  write.table(table, paste("humann_keggAsCol_withRurUrb_allkegg_", names[i], ".txt", sep=""), 
              sep="\t", row.names = F, col.names = T, quote=F)
  print(ncol(table)) #223 for mpm, 287 for mpt
  
  ##remove all columns that are all 0
  table = table[,c(T, T, colSums(table[,-(1:2)]) > 0)]
  print(ncol(table)) #215 for mpm, 287 for mpt
  write.table(table, paste("humann_keggAsCol_withRurUrb_", names[i], ".txt", sep=""), 
              sep="\t", row.names = F, col.names = T, quote=F)
  
  ##log
  count = 0.00001
  for(c in 7:ncol(table)) {
    for(r in 1:nrow(table)) {
      table[r,c] = log10(table[r,c] + count) + -log10(count)
    }
  }
  write.table(table, paste("humann_keggAsCol_log_", names[i], ".txt", sep=""), 
              sep="\t", row.names = F, col.names = T, quote=F)
}