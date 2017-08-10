##rarefy 16S data, then recalculate diversity measures
##plot the average of 10 repetitions

rm(list=ls())
setwd("")

library(vegan)
library(nlme)

levels = c("phylum","class","order","family","genus", "otu")
numReps = 10

##get metadata
meta = read.table("phylum_taxaAsColumnsLogNorm_withMetadata.txt", header=T, sep="\t", colClasses = "character")
meta = meta[,1:5]

##caclulate the p-values using mixed models based on the given sampleIDs and diversity (div)
##write results to the given row in the given p table, which is returned
calculateP <- function(sampleID, div, pTab, row) {
  df = data.frame(sampleID, div)
  table = merge(meta, df, by="sampleID")
  measure = table$div
  time = factor(table$timepoint)
  urbanRural = factor(table$ruralUrban)
  patientID = factor(table$patientID)
  
  reducedModel <- gls(measure ~ time + urbanRural , method="REML")
  
  fullModelLME <- lme(measure ~ time + urbanRural , method="REML", random = ~1|patientID)		
  
  pTab$pTime[row] <- anova(fullModelLME)$"p-value"[2]
  pTab$pUrbanRural[row] <- anova(fullModelLME)$"p-value"[3]
  pTab$pSubject[row] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
  return(pTab)
}

##set up p-value tables
pRich = data.frame(level = levels, 
                     pUrbanRural = rep(NA, length(levels)),
                     pTime = rep(NA, length(levels)),
                     pSubject = rep(NA, length(levels)),
                     stringsAsFactors = F)
pShan = data.frame(level = levels, 
                   pUrbanRural = rep(NA, length(levels)),
                   pTime = rep(NA, length(levels)),
                   pSubject = rep(NA, length(levels)),
                   stringsAsFactors = F)
pInvSimp = data.frame(level = levels, 
                   pUrbanRural = rep(NA, length(levels)),
                   pTime = rep(NA, length(levels)),
                   pSubject = rep(NA, length(levels)),
                   stringsAsFactors = F)
pEven = data.frame(level = levels, 
                   pUrbanRural = rep(NA, length(levels)),
                   pTime = rep(NA, length(levels)),
                   pSubject = rep(NA, length(levels)),
                   stringsAsFactors = F)

jpeg("rarefyDiversity.jpg", height=4200, width=3650, res=300)
layout(matrix(1:24, nrow=4, ncol=6, byrow=F), heights=1, widths=1)
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.9, oma=c(2,0,0,0), mar=c(4.1, 6, 4,.1))

##run
for(lev in levels) {
  print(lev)
  ##read table
  file = paste(lev, "_taxaAsColumns.txt", sep="")
  if(lev == "otu") {
    file = "abundantOTUForwardTaxaAsColumns.txt"
  }
  table = read.table(file, header=T, sep="\t")
  
  ##fix sampleID
  table$sample = as.character(table$sample)
  names(table)[1] = "sampleID"
  
  ##fix OTU sampleID
  if(lev == "otu") {
    table$sampleID = paste(table$sampleID, "_1", sep="")
  }
  
  ##forward reads only
  table = table[grepl("_1", table$sampleID),]
  
  ##use sample with smallest number of reads as subsample size
  numReads = rowSums(table[,-1])
  sub = min(numReads)
  print(sub)
  
  ##10 repetitions of rarefecation + diversity measures
  aveRich = rep(0, nrow(table))
  aveShan = rep(0, nrow(table))
  aveInvSimp = rep(0, nrow(table))
  aveEven = rep(0, nrow(table))
  for(i in 1:numReps) {
    ##rarefy
    rar = rrarefy(table[,-1], sub)
    rich = specnumber(rar)
    shan = diversity(rar, index="shannon")
    invSimp = diversity(rar, index="invsimpson")
    even = shan / log(rich)
    
    ##write results
    write.table(data.frame(sampleID = table$sample, richness=rich, shannon=shan, invSimpson=invSimp, evenness=even),
                paste("rarefyDiversity_diversityValues_", lev, "_rep", i, ".txt", sep=""),
                sep="\t", row.names = F, col.names = T, quote=F)
    write.table(cbind(sampleID=table$sample, rar),
                paste("rarefyDiversity_rarefiedCounts_", lev, "_rep", i, ".txt", sep=""),
                sep="\t", row.names = F, col.names = T, quote=F)
    
    ##add to averages
    aveRich = aveRich + rich 
    aveShan = aveShan + shan
    aveInvSimp = aveInvSimp + invSimp
    aveEven = aveEven + even
  }
  
  ##get average
  aveRich = aveRich / numReps
  aveShan = aveShan / numReps
  aveInvSimp = aveInvSimp / numReps
  aveEven = aveEven / numReps
  results = data.frame(sampleID = table$sampleID, richness=aveRich, shannon=aveShan, 
                       invSimpson=aveInvSimp, evenness=aveEven)
  write.table(results,
              paste("rarefyDiversity_diversityValues_", lev, "_average.txt", sep=""),
              sep="\t", row.names = F, col.names = T, quote=F)
  
  ##p-values for average
  row = lev == levels
  sID = results$sampleID
  pRich = calculateP(sID, aveRich, pRich, row)
  pShan = calculateP(sID, aveShan, pShan, row)
  pInvSimp = calculateP(sID, aveInvSimp, pInvSimp, row)
  pEven = calculateP(sID, aveEven, pEven, row)
  
  ##plot average
  mrg = merge(meta, results, by="sampleID")
  urb = factor(mrg$ruralUrban)
  x=ifelse(urb=="rural", 1, 2)
  col=ifelse(urb=="rural", "blue", "red")
  points=ifelse(mrg$timepoint=="first_A", 16, 17)
  
  boxplot(mrg$richness~urb, ylab="Richness") 
          # main=paste(lev, "\np=", format(pRich$pUrbanRural[row], digits = 3), sep=""))
  points(mrg$richness~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pRich$pUrbanRural[row], digits=2),sep=""), side=3, line=0)
  mtext(lev, side = 3, line = 2, cex=1.5)
  
  boxplot(mrg$shannon~urb, ylab="Shannon Index")
  points(mrg$shannon~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pShan$pUrbanRural[row], digits=2),sep=""), side=3, line=0)
  
  boxplot(mrg$invSimpson~urb, ylab="Inverse Simpson Index")
  points(mrg$invSimpson~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pInvSimp$pUrbanRural[row], digits=2),sep=""), side=3, line=0)
  
  boxplot(mrg$evenness~urb, ylab="Evenness")
  points(mrg$evennes~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pEven$pUrbanRural[row], digits=2),sep=""), side=3, line=0)
}

##legend
par(oma=c(0,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", cex=1.5, horiz=T,
       legend=c("rural", "urban"),
       col=c("blue", "red"),
       pch=16)
dev.off()

##adjust p-values
##function that corrects all p-values for the given table
correctP <- function(pTab) {
  pTab$adjustedPurbanRural = p.adjust(pTab$pUrbanRural, method="BH")
  pTab$adjustedPtime = p.adjust(pTab$pTime, method="BH")
  pTab$adjustedPsubject = p.adjust(pTab$pSubject, method="BH")
  return(pTab)
}
pRich = correctP(pRich)
pShan = correctP(pShan)
pInvSimp = correctP(pInvSimp)
pEven = correctP(pEven)

##write p-values
write.table(pRich, "rarefyDiversity_pValues_richness.txt", sep="\t", row.names = F, col.names = T, quote=F)
write.table(pShan, "rarefyDiversity_pValues_shannon.txt", sep="\t", row.names = F, col.names = T, quote=F)
write.table(pInvSimp, "rarefyDiversity_pValues_invsimpson.txt", sep="\t", row.names = F, col.names = T, quote=F)
write.table(pEven, "rarefyDiversity_pValues_evenness.txt", sep="\t", row.names = F, col.names = T, quote=F)