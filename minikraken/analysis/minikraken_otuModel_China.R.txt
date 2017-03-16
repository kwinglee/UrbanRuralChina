##models for minikraken taxa

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")

setwd("")

taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")

for(taxa in taxaLevels) {
  print(taxa)
  file = paste("minikraken_merged_taxaAsCol_logNorm_", taxa, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  names <- vector()
  pValuesUrbanRural <- vector()
  meanBug <- vector()
  sdBug <- vector()
  meanUrban <- vector()
  sdUrban <- vector()
  meanRural <- vector()
  sdRural <- vector()
  pValuesUrbanRuralWilcox <- vector()
  r.squared <- vector()
  index <- 1
  
  ##p-values
  for(i in 4:ncol) {
    if(sum(table[,i] != 0 ) > nrow(table) / 4) {
      
      bug <- table[,i]
      meanBug[index] <- mean(bug)
      meanUrban[index] <- mean(bug[table$ruralUrban=="urban"])
      meanRural[index] <- mean(bug[table$ruralUrban=="rural"])
      sdBug[index] <- sd(bug)
      sdUrban[index] <- sd(bug[table$ruralUrban=="urban"])
      sdRural[index] <- sd(bug[table$ruralUrban=="rural"])
      urbanRural <- factor(table$ruralUrban)
      
      ##linear model
      model = lm(bug~urbanRural)
      pValuesUrbanRural[index] <- anova(model)$`Pr(>F)`[1]
      r.squared[index] = summary(model)$r.squared
      names[index] = names(table)[i]
      
      ##non parametric test
      pValuesUrbanRuralWilcox[index] = wilcox.test(bug~urbanRural)$p.value
      
      index=index+1
      
    }
  }
  dFrame <- data.frame(names, meanBug, meanUrban, sdUrban, meanRural, sdRural,
                       pValuesUrbanRural, pValuesUrbanRuralWilcox)
  dFrame$UrbanToRural <- meanUrban / meanRural
  dFrame$adjustedPurbanRural <- p.adjust( dFrame$pValuesUrbanRural, method = "BH" )
  dFrame$adjustedPurbanRuralWilcox <- p.adjust(dFrame$pValuesUrbanRuralWilcox, method="BH")
  dFrame$r.squared = r.squared
  dFrame <- dFrame [order(dFrame$pValuesUrbanRural),]
  write.table(dFrame, file=paste("minikraken_otuModel_pValues_", taxa, ".txt", sep=""), sep="\t",row.names=FALSE)
  
  ##plots
  pdf(paste("kraken_model_boxplots_", taxa, ".pdf", sep=""), height=5, width=5)
  urbanRural = factor(table$ruralUrban)
  for(i in 1:nrow(dFrame)) {
    name = dFrame$names[i]
    abun = table[,names(table)==name]
    
    graphMain =  paste("WGS ", taxa, ": ", name,
                       "\npAdjRuralUrban= ", format(dFrame$adjustedPurbanRural[i],digits=3), sep="")
    boxplot(abun~urbanRural, main=graphMain, ylab="log relative abundance", cex.main=1, outline=F, ylim=range(abun))
    points(abun~jitter(as.numeric(urbanRural)), pch=16, col=ifelse(urbanRural=="rural", "blue", "red"))
  }
  dev.off()
}