##linear model on metabolon and metadata pcoa axes

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")

setwd("")

anal = c("metabolon", "metadata")

for(a in anal) {
  fileName = paste("pcoaCorrected_", a, ".txt", sep="")
  table <-read.table(fileName,header=TRUE,sep="\t")
  numCols <- ncol(table)
  colClasses <- c("character", "character", rep("numeric", numCols-2))
  table <-read.table(fileName,header=TRUE,sep="\t",colClasses=colClasses)
              
  ##output vectors
  names <- vector()
  pValuesRuralUrban <- vector()
  meanMDS <- vector()
  meanUrban <- vector()
  meanRural <- vector()
  index <- 1
  pdf(paste("pcoaModel_boxplots_", a, ".pdf", sep=""))
  
  for( i in 3:ncol(table)) {
    if(sum(table[,i] != 0 ) > nrow(table) / 4  && !all(table[,i]==1)) { #only analyze if more than a quarter are nonzero
      #in metabolon table three columns (173, 260, 277) give a gls error: computed "gls" fit is singular, rank 2 -> these columns are all 1
      mds <- table[,i]
      meanMDS[index] <- mean(mds)
      meanUrban[index] <- mean(mds[table$ruralUrban=="urban"])
      meanRural[index] <- mean(mds[table$ruralUrban=="rural"])
      ruralUrban <- factor(table$ruralUrban)
      
      myFrame <- data.frame(mds, ruralUrban )
      
      model <- gls(mds ~ ruralUrban, method="REML",	data = myFrame )
      
      pValuesRuralUrban[index] <- anova(model)$"p-value"[2]
      names[index] = names(table)[i]
      
      graphMain =  paste(names(table)[i], "\n",
                         " pRuralUrban= ", format(pValuesRuralUrban[index],digits=3), sep="")
      boxplot(mds ~ ruralUrban, main=graphMain, xlab="", ylab="scaled value")
      points(x=ruralUrban, y=mds, col=ifelse(table$ruralUrban == "rural", "blue", "red"), pch=16)
      
      index=index+1
      
    }
  }
  
  dFrame <- data.frame( names,meanMDS,meanUrban, meanRural, pValuesRuralUrban)
  dFrame$UrbanToRural <- meanUrban/meanRural
  dFrame$adjustedPruralUrban <- p.adjust( dFrame$pValuesRuralUrban, method = "BH" )
  write.table(dFrame, file=paste("pcoaModel_pValues_", a, ".txt", sep=""), sep="\t",row.names=FALSE, quote=F)
  dev.off()
}