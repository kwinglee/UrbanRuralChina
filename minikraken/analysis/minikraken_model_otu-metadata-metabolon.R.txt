##find correlations between minikraken, metadata and metabolon data

rm(list=ls())
library(Kendall)

setwd("")

metadata = read.table(paste("metadata_cleaned11-2015.txt", sep=""), 
                      header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 55)))
metabolon = read.table(paste("MetabolonScaledImpData-transpose.txt", sep=""), 
                       header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 337)))

categoricalVars = c("F9", "F13", "F17", "sex") #categorical metadata variables

##get sample ids to be the same -> convert to numeric
names(metabolon)[1] = "sampleID"
names(metabolon)[2] = "ruralUrban"

##function that runs the comparisons
##table = table to compare, splitCol = split point between datasets (last column in OTU dataset)
##name = name of output file
runComparisons <- function(table, splitCol, name) {
  namesOTU <- vector()
  namesOther <- vector()
  r <- vector()
  r.squared <- vector()
  pValues <- vector()
  pNonparam <- vector()
  nonparamTest <- vector()
  index <- 1
  
  pdf(paste("minikraken_model_plots_", name, ".pdf",sep=""))
  
  for(i in 4:splitCol) {
    if(sum( table[,i] > 0) > nrow(table) /4 && !all(table[,i]==1)) {
      for ( j in (splitCol+1):ncol(table)) {
        if(!all(table[,j]==1)) {
          namesOTU[index] <- names(table)[i]
          namesOther[index] <- names(table)[j]
          
          if(namesOTU[index] %in% categoricalVars | namesOther[index] %in% categoricalVars) {
            valI = table[,i]
            valJ = table[,j]
            if(namesOTU[index] %in% categoricalVars) {
              valI = table[,j]
              valJ = table[,i]
            }
            valJ = factor(valJ)
            
            r[index] = NA
            pValues[index] = t.test(valI~valJ)$p.value
            r.squared[index] = summary(lm(valI~valJ))$r.squared
            
            pNonparam[index] = wilcox.test(valI~valJ)$p.value
            nonparamTest[index] = "Wilcox"
            
            if(!is.na(pNonparam[index]) & pNonparam[index] < 0.05) { 
              myText <- paste(namesOTU[index], " vs. " ,namesOther[index] ,"\n", "p=" ,  
                              format(pValues[index] , digits=3), "r=", 
                              format(r[index], digits=3), "kendall p=" , 
                              format(pNonparam[index], digits=3)) 
              boxplot(valI~valJ, main=myText, xlab=namesOther[index], ylab=namesOTU[index])
              points(valI~valJ, pch=16, col = ifelse(table$ruralUrban == "rural", "blue", "red"))
            }
            
          } else {
            r[index] <- cor(table[,i], table[,j], method="spearman")
            aLm <- lm(table[,i] ~ table[,j])
            pValues[index] <- anova(aLm)$"Pr(>F)"[1]
            r.squared[index] = summary(aLm)$r.squared
            pNonparam[index] <- Kendall(table[,i], table[,j])$sl[1]
            nonparamTest[index] <- "Kendall"
            
            if((!is.na(pNonparam[index]) & pNonparam[index] < 0.05) |
               (pValues[index] < 0.05)) { 
              myText <- paste(namesOTU[index], " vs. " ,namesOther[index] ,"\n", "p=" ,  
                              format(pValues[index] , digits=3), "r=", 
                              format(r[index], digits=3), "kendall p=" , 
                              format(pNonparam[index], digits=3)) 
              plot(table[,j],table[,i] , main=myText, pch=16, xlab=namesOther[index], ylab=namesOTU[index],
                   col = ifelse(table$ruralUrban == "rural", "blue", "red"))
              abline(aLm)
            }
          }
          
          index <- index + 1
        }
      }
    }
  }
  
  dev.off()
  dFrame <- data.frame( namesOTU, namesOther, pLM=pValues, r, r.squared, nonparamTest, pNonparam) 
  dFrame <- dFrame [order(dFrame$pLM),]
  dFrame$adjLMp <-  p.adjust(dFrame$pLM, method = "BH")
  dFrame$adjNonparmP <- p.adjust(dFrame$pNonparam, method = "BH")
  
  write.table( file= paste("minikraken_model_pValues_", name, ".txt", sep=""), dFrame, row.names=FALSE, sep="\t")
}


taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")

for(taxa in taxaLevels ) {
  print(taxa)
  file = paste("minikraken_merged_taxaAsCol_logNorm_", taxa, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  ##convert sample IDs to numeric like metadata and metabolon
  names = table$sampleID
  names = sub("A", "", names)
  table$sampleID = as.numeric(names)
  
  ##metadata vs OTU
  print("  metadata")
  mrg.md = merge(table, metadata, by=c("sampleID", "ruralUrban"))
  runComparisons(mrg.md, ncol, paste("metadata_v_", taxa, sep=""))
  
  ##metabolite dataset
  print("  metabolite")
  mrg.mb = merge(table, metabolon, by=c("sampleID", "ruralUrban"))
  runComparisons(mrg.mb, ncol, paste("metabolon_v_", taxa, sep=""))
}
