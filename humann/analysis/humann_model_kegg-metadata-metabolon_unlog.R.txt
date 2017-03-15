##find correlations between unlogged HUMAnN, metadata and metabolon data

rm(list=ls())
library("Kendall")

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
##name = name of output file, desc = table with descriptions of the KEGG data
runComparisons <- function(table, splitCol, name, desc) {
  namesKEGG <- vector()
  description <- vector()
  namesOther <- vector()
  r <- vector()
  r.squared <- vector()
  pValues <- vector()
  pNonparam <- vector()
  nonparamTest <- vector()
  index <- 1
  
  pdf(paste("humann_model_unlog_", name, "_plots.pdf",sep=""))
  
  for(i in 7:splitCol) {
    if(sum(table[,i] >0 ) > nrow(table) /4 && length(unique(table[,i])) > 1) {
      for ( j in (splitCol+1):ncol(table)) {
        if(!all(table[,j]==1)) { 
          namesKEGG[index] <- names(table)[i]
          namesOther[index] <- names(table)[j]
          description[index] <- desc$NAME[desc$sampleID==namesKEGG[index]]
          
          if(namesKEGG[index] %in% categoricalVars | namesOther[index] %in% categoricalVars) {
            valI = table[,i]
            valJ = table[,j]
            if(namesKEGG[index] %in% categoricalVars) {
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
              myText <- paste(namesKEGG[index], " vs. " ,namesOther[index] ,"\n", "p=" ,  
                              format(pValues[index] , digits=3), "r=", 
                              format(r[index], digits=3), "kendall p=" , 
                              format(pNonparam[index], digits=3)) 
              boxplot(valI~valJ, main=myText, xlab=namesOther[index], ylab=namesKEGG[index])
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
              myText <- paste(namesKEGG[index], " vs. " ,namesOther[index] ,"\n", "p=" ,  
                              format(pValues[index] , digits=3), "r=", 
                              format(r[index], digits=3), "kendall p=" , 
                              format(pNonparam[index], digits=3)) 
              plot(table[,j],table[,i] , main=myText, pch=16, xlab=namesOther[index], ylab=namesKEGG[index],
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
  dFrame <- data.frame(namesKEGG, description, namesOther, pLM=pValues, 
                       r, r.squared, nonparamTest, pNonparam) 
  dFrame <- dFrame [order(dFrame$pLM),]
  dFrame$adjLMp <- p.adjust(dFrame$pLM, method = "BH")
  dFrame$adjNonparmP <- p.adjust(dFrame$pNonparam, method = "BH")
  
  write.table(file= paste("humann_model_unlog_", name, "_pValues.txt", sep=""), dFrame, 
               row.names=FALSE, sep="\t", quote = F)
}

levels = c("module", "pathway")

for(lev in levels) {
  print(lev)
  file = paste("humann_keggAsCol_withRurUrb_", lev, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", ncol-2)))
  
  desc = read.table(paste("humann_keggAsRow_", lev, ".txt", sep=""), sep="\t", quote="", 
                    header=T, stringsAsFactors = F)
  
  ##convert sample IDs to numeric like metadata and metabolon
  names = table$sampleID
  names = sub("A", "", names)
  table$sampleID = as.numeric(names)
  
  ##metadata vs OTU
  mrg.md = merge(table, metadata, by=c("sampleID", "ruralUrban"))
  print(" metadata")
  runComparisons(mrg.md, ncol, paste("metadata_v_", lev, sep=""), desc)
  
  ##metabolite dataset
  mrg.mb = merge(table, metabolon, by=c("sampleID", "ruralUrban"))
  print(" metabolite")
  runComparisons(mrg.mb, ncol, paste("metabolon_v_", lev, sep=""), desc)
}
