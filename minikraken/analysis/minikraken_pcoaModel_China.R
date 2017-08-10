##models for minikraken PCoA axes

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")


setwd("")

taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")

for(taxa in taxaLevels) 
{
  print(taxa)
  file = paste("minikraken_PcoaCorrected_", taxa, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  names <- vector()
  pValuesUrbanRural <- vector()
  meanBug <- vector()
  meanUrban <- vector()
  meanRural <- vector()
  pValuesUrbanRuralWilcox <- vector()
  index <- 1
  pdf(paste("minikraken_pcoaModel_boxplots_", taxa, ".pdf", sep=""))
  
  for(i in 4:ncol)
    if(sum(table[,i] != 0 ) > nrow(table) / 4 )
    {
      
      bug <- table[,i]
      meanBug[index] <- mean(bug)
      meanUrban[index] <- mean(bug[table$ruralUrban=="urban"])
      meanRural[index] <- mean(bug[table$ruralUrban=="rural"])
      urbanRural <- factor(table$ruralUrban)
      
      ##linear model
      pValuesUrbanRural[index] <- anova(lm(bug~urbanRural))$`Pr(>F)`[1]
      names[index] = names(table)[i]
      
      ##non parametric test
      pValuesUrbanRuralWilcox[index] = wilcox.test(bug~urbanRural, exact=F)$p.value
      #note, if don't give exact=F, get warning: cannot compute exact p-value with ties
      
      ##plot
      graphMain =  paste(names[index],"\npRuralUrban= ", format(pValuesUrbanRural[index],digits=3), sep="")
      boxplot(bug~urbanRural, main=graphMain, ylab="MDS", cex.main=1)
      points(bug~urbanRural, pch=16, col=ifelse(urbanRural=="rural", "blue", "red"))
      
      index=index+1
      
    }
  
  dFrame <- data.frame(names, meanBug, meanUrban, meanRural, pValuesUrbanRural, pValuesUrbanRuralWilcox)
  dFrame$UrbanToRural <- meanUrban / meanRural
  dFrame$adjustedPurbanRural <- p.adjust( dFrame$pValuesUrbanRural, method = "BH" )
  dFrame$adjustedPurbanRuralWilcox <- p.adjust(dFrame$pValuesUrbanRuralWilcox, method="BH")
  write.table(dFrame, file=paste("minikraken_pcoaModel_pValues_", taxa, ".txt", sep=""), sep="\t",row.names=FALSE)
  dev.off()
}