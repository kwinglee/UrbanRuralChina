##models for HUMAnN PCoA

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")

setwd("")

levels = c("module", "pathway")

for(lev in levels ) 
{
  print(lev)
  file=paste("humann_PcoaCorrected_unlog_", lev, ".txt",sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", ncol-2)))
  
  names <- vector()
  pValuesUrbanRural <- vector()
  meanPCoA <- vector()
  meanUrban <- vector()
  meanRural <- vector()
  pValuesUrbanRuralWilcox <- vector()
  index <- 1
  pdf(paste("humann_pcoaModel_boxplots_unlog_", lev, ".pdf", sep=""))
  
  for(i in 3:ncol)
    if(sum(table[,i] != 0 ) > nrow(table) / 4 )
    {
      
      pc <- table[,i]
      meanPCoA[index] <- mean(pc)
      meanUrban[index] <- mean(pc[table$ruralUrban=="urban"])
      meanRural[index] <- mean(pc[table$ruralUrban=="rural"])
      urbanRural <- factor(table$ruralUrban)
      
      ##linear model
      pValuesUrbanRural[index] <- anova(lm(pc~urbanRural))$`Pr(>F)`[1]
      names[index] = names(table)[i]
      
      ##non parametric test
      pValuesUrbanRuralWilcox[index] = wilcox.test(pc~urbanRural, exact=F)$p.value
      #note, if don't give exact=F, get warning: cannot compute exact p-value with ties
      
      ##plot
      graphMain =  paste(names[index],"\npRuralUrban= ", format(pValuesUrbanRural[index],digits=3), sep="")
      boxplot(pc~urbanRural, main=graphMain, ylab="MDS", cex.main=.6)
      points(pc~urbanRural, pch=16, col=ifelse(urbanRural=="rural", "blue", "red"))
      
      index=index+1
      
    }
  
  dFrame <- data.frame(names, meanPCoA, meanUrban, meanRural, pValuesUrbanRural, pValuesUrbanRuralWilcox)
  dFrame$UrbanToRural <- meanUrban / meanRural
  dFrame$adjustedPurbanRural <- p.adjust( dFrame$pValuesUrbanRural, method = "BH" )
  dFrame$adjustedPurbanRuralWilcox <- p.adjust(dFrame$pValuesUrbanRuralWilcox, method="BH")
  write.table(dFrame, file=paste("humann_pcoaModel_pValues_unlog_", lev, ".txt",sep=""), sep="\t",row.names=FALSE)
  dev.off()
}