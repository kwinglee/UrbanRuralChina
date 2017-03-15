##linear model between 16S, WGS, metabolite and metadata to proportion of reads mapped to virulence databases

rm(list=ls())
setwd("")
library(Kendall)

databases = c("VFDBcore", "VFDBfull", "MvirDB")
ncols = c(1059, 9047, 15263)
start = 5

##for the given table, generates associations between the proportion of reads that aligned to the virulence database
##and the table
##name is the suffix to use, if desc is not NA (humann), include a description, vir is the database, db is the database name
assocPropReads <- function(table, name, vir, db, desc=NA) {
  nvir = ncol(vir)
  mrg = merge(vir, table, by="sampleID")
  prop = mrg$proportionReadsMapped
  
  namesOther = vector()
  r = vector()
  r.squared = vector()
  pValues = vector()
  kendallP = vector()
  description = vector()
  index = 1
  
  pdf(paste("virulenceAssoc_plots_", db, "_v_", name, ".pdf", sep=""))
  for(i in (nvir+1):ncol(mrg)) {
    namesOther[index] = names(mrg)[i]
    if(!is.na(desc)) {
      description[index] <- desc$NAME[desc$sampleID==namesOther[index]]
    }
    
    r[index] = cor(prop, mrg[,i], method="spearman")
    aLm = lm(prop ~ mrg[,i])
    pValues[index] = anova(aLm)$"Pr(>F)"[1]
    r.squared[index] = summary(aLm)$r.squared
    kendallP[index] = Kendall(prop, mrg[,i])$sl[1]
    
    myText = paste("Proportion of reads mapped vs. " ,namesOther[index] ,"\n", "p=" ,  
                   format(pValues[index] , digits=3), "r=", 
                   format(r[index], digits=3), "kendall p=" , 
                   format(kendallP[index], digits=3)) 
    plot(mrg[,i], prop, main=myText, pch=16, xlab=namesOther[index], ylab="Proportion of reads mapped",
         col = ifelse(mrg$ruralUrban == "rural", "blue", "red"))
    abline(aLm)
    
    index = index + 1
  }
  dev.off()
  
  if(is.na(desc)) {
    dFrame = data.frame(namesOther, kendallP, pLM=pValues, r, r.squared) 
  } else {
    dFrame = data.frame(namesOther, description, kendallP, pLM=pValues, r, r.squared) 
  }
  dFrame = dFrame[order(dFrame$pLM),]
  dFrame$adjKendallP =  p.adjust(dFrame$kendallP, method = "BH")
  dFrame$adjLMp =  p.adjust(dFrame$pLM, method = "BH")
  
  write.table(file= paste("virulenceAssoc_pValues_", db, "_v_", name, ".txt", sep=""), dFrame, row.names=FALSE, sep="\t")
  
}

##for the given table, generates associations between the proportion of reads that aligned to each virulence gene
##and the table
##name is the suffix to use, if desc is not NA (humann), include a description, db is the name of the database
assocGenes <- function(table, name, vir, db, desc=NA) {
  nvir = ncol(vir)
  mrg = merge(vir, table, by="sampleID")
  
  length = (nvir-2-start) * (ncol(mrg)-nvir)
  namesVirulence = rep(NA, length)
  namesOther = rep(NA, length)
  r = rep(NA, length)
  r.squared = rep(NA, length)
  pValues = rep(NA, length)
  kendallP = rep(NA, length)
  if(!is.na(desc)) {
    description = rep(NA, length)
  }
  
  index = 1
  
  pdf(paste("virulenceAssoc_allGenes_plots_", db, "_v_", name, ".pdf",sep=""))
  
  for( i in (start+2):nvir) {
    if(!grepl("human", names(mrg)[i]) && !grepl("Human", names(mrg)[i]))  { #don't use human filtered
      for ( j in (nvir+1):ncol(mrg)) {
        namesVirulence[index] = names(mrg)[i]
        namesOther[index] = names(mrg)[j]
        if(!is.na(desc)) {
          description[index] <- desc$NAME[desc$sampleID==namesOther[index]]
        }
        
        r[index] = cor(mrg[,i], mrg[,j], method="spearman")
        aLm = lm(mrg[,i] ~ mrg[,j])
        r.squared[index] = summary(aLm)$r.squared
        pValues[index] = anova(aLm)$"Pr(>F)"[1]
        kendallP[index] = Kendall(mrg[,i], mrg[,j])$sl[1]
        
        if(pValues[index] <= 0.05) { #reduce size of pdfs by only plotting potentially significant associations
          myText = paste(namesVirulence[index], " vs. ", namesOther[index] ,"\n", "p=" ,  
                         format(pValues[index] , digits=3), "r=", 
                         format(r[index], digits=3), "kendall p=" , 
                         format(kendallP[index], digits=3)) 
          plot(mrg[,j],mrg[,i] , main=myText, pch=16, xlab=namesOther[index], ylab=namesVirulence[index],
               col = ifelse(mrg$ruralUrban == "rural", "blue", "red"))
          abline(aLm)
        }
        
        index = index + 1
      }
    }
  }
  
  dev.off()
  if(!is.na(desc)) {
    dFrame = data.frame(namesVirulence, description, namesOther, kendallP, pLM=pValues, r, r.squared) 
  } else {
    dFrame = data.frame(namesVirulence, namesOther, kendallP, pLM=pValues, r, r.squared) 
  }
  dFrame = dFrame[1:(index-1),]
  dFrame = dFrame [order(dFrame$pLM),]
  dFrame$adjKendallP =  p.adjust(dFrame$kendallP, method = "BH")
  dFrame$adjLMp=  p.adjust(dFrame$pLM, method = "BH")
  
  write.table( file= paste("virulenceAssoc_allGenes_pValues_", db, "_v_", name, ".txt", sep=""), dFrame, row.names=FALSE, sep="\t")
}

for(d in 1:length(databases)) {
  db = databases[d]
  print(db)
  fname = paste(db, "_results.txt", sep="")
  vir = read.table(fname, header=T, sep="\t", quote="", colClasses=c("character", rep("numeric", ncols[d]-1)))
  
  names(vir)[1] = "sampleID"
  filter <- function(i) (return(sum(vir[,i] != 0 ) > nrow(vir) / 4))
  keep = sapply(1:ncol(vir), filter)
  vir = vir[,keep]
  
  ##add metadata
  meta = read.table("genus_taxaAsColumns_relAbunFwd.txt", sep="\t", header=T,
                    colClasses=c(rep("character", 3), rep("numeric", 347)))
  meta = meta[meta$timepoint=="first_A",1:2]
  meta$sampleID = gsub("_1", "", meta$sampleID)
  vir = merge(meta, vir, by="sampleID")
  nvir = ncol(vir)
  
  #####16S
  print(" 16S")
  ##get OTU taxonomy for the given OTU name
  taxonomy = read.table("abundantOTU.chinaForward.taxonomy.txt", header=F, colClasses="character", sep="\t")
  ##function that returns full taxonomy for the given name in the given level
  getOTUTaxonomy <- function(name) {
    n = sub("X", "", name)
    n = paste("Consensus", n, sep="")
    r = taxonomy$V1==n
    return(paste("p_", taxonomy$V6[r], ";c_", taxonomy$V9[r], ";o_", taxonomy$V12[r], ";f_", taxonomy$V15[r], ";g_", taxonomy$V18[r], ";", n, sep=""))
  }
  
  taxaLevels = c("phylum","class","order","family","genus", "otu")
  for(taxa in taxaLevels) {
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
    table = table[table$readNumber == "1",]
    
    ##second time point will be removed when merge table
    ##fix sample IDs
    table$sampleID = sub("_1", "", table$sampleID)
    
    ##remove metadata
    table = table[,-(2:5)]
    
    ##filter zeros in > 1/4 samples
    filter <- function(i) (return(sum(table[,i] != 0 ) > nrow(table) / 4))
    keep = sapply(1:ncol(table), filter)
    table = table[,keep]
    
    ##get description for OTU
    if(taxa=="otu") {
      id = names(table)[-1]
      names = sapply(id, getOTUTaxonomy)
      desc = data.frame(sampleID=id, NAME=names, stringsAsFactors = F)
    }
    
    ##association to proportion of reads
    name = paste("16S_", taxa, sep="")
    if(taxa=="otu") {
      assocPropReads(table, name, vir, db, desc)
    } else {
      assocPropReads(table, name, vir, db)
    }
    
    ##association to individual genes
    if(taxa=="otu") {
      assocGenes(table, name, vir, db, desc)
    } else {
      assocGenes(table, name, vir, db)
    }
  }
  
  #####minikraken
  print(" minikraken")
  taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")
  for(taxa in taxaLevels) {
    print(taxa)
    file = paste("minikraken_merged_taxaAsCol_logNorm_", taxa, ".txt", sep="")
    table = read.table(file, header=T, sep="\t")
    ncol = ncol(table)
    table = read.table(file, header=T, sep="\t", 
                       colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
    
    ##remove metadata
    table = table[,-(2:3)]
    
    ##filter zeros in > 1/4 samples
    filter <- function(i) (return(sum(table[,i] != 0 ) > nrow(table) / 4))
    keep = sapply(1:ncol(table), filter)
    table = table[,keep]
    
    ##association to proportion of reads
    name = paste("minikraken_", taxa, sep="")
    assocPropReads(table, name, vir, db)
    
    ##association to individual genes
    assocGenes(table, name, vir, db)
  }
  
  #####humann
  print(" humann")
  levels = c("module", "pathway")
  for(lev in levels) {
    print(lev)
    file = paste("humann_keggAsCol_withRurUrb_", lev, ".txt", sep="")
    table = read.table(file, header=T, sep="\t")
    ncol = ncol(table)
    table = read.table(file, header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", ncol-2)))
    
    desc = read.table(paste("humann_keggAsRow_", lev, ".txt", sep=""), sep="\t", quote="", 
                      header=T, stringsAsFactors = F)
    
    ##remove metadata
    table = table[,-(2:6)]
    
    ##filter zeros in > 1/4 samples
    filter <- function(i) (return(sum(table[,i] != 0 ) > nrow(table) / 4))
    keep = sapply(1:ncol(table), filter)
    table = table[,keep]
    
    ##association to proportion of reads
    name = paste("humann_", lev, sep="")
    assocPropReads(table, name, vir, db, desc)
    
    ##association to individual genes
    assocGenes(table, name, vir, db, desc)
  }
  
  #####metadata
  print(" metadata")
  metadata = read.table("metadata_cleaned11-2015.txt", 
                        header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 55)))
  ##remove ruralUrban
  metadata = metadata[,-2]
  ##fix sampleIDs
  metadata$sampleID[metadata$sampleID > 80 & metadata$sampleID < 100] = paste("0", metadata$sampleID[metadata$sampleID > 80 & metadata$sampleID < 100], sep="")
  metadata$sampleID = paste(metadata$sampleID, "A", sep="")
  
  name = "metadata"
  assocPropReads(metadata, name, vir, db)
  assocGenes(metadata, name, vir, db)
  
  #####metabolon
  print(" metabolon")
  metabolon = read.table("MetabolonScaledImpData-transpose.txt", 
                         header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 337)))
  
  ##remove ruralUrban
  metabolon = metabolon[,-2]
  ##fix sampleIDs
  names(metabolon)[1] = "sampleID"
  metabolon$sampleID[metabolon$sampleID > 80 & metabolon$sampleID < 100] = 
    paste("0", metabolon$sampleID[metabolon$sampleID > 80 & metabolon$sampleID < 100], sep="")
  metabolon$sampleID = paste(metabolon$sampleID, "A", sep="")
  ##remove the metabolites that are all 1
  filter <- function(i) (return(!all(metabolon[,i] == 1)))
  keep = sapply(1:ncol(metabolon), filter)
  metabolon = metabolon[,keep]
  
  name = "metabolon"
  assocPropReads(metabolon, name, vir, db)
  assocGenes(metabolon, name, vir, db)
}