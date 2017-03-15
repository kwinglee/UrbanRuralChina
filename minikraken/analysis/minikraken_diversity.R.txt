##analysis of diversity measures for minikraken results

rm(list=ls())

library(vegan)

setwd("")

taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")

##draws the plots for richness, inverse simpson, shannon and shannon evenness for the given table
##uses the given taxa name and writes the measures to the file outputName, puts the values in pVal
drawPlots <- function(table, taxa, outputName, pVal) {
  urb = factor(table$ruralUrban)
  col=ifelse(urb=="rural", "blue", "red")
  points=16
  x=ifelse(urb=="rural", 1, 2)
  row = pVal$level==taxa
  
  ##richness
  rich = specnumber(table[,-(1:3)])
  boxplot(rich~urb, ylab="Richness")
  points(rich~jitter(x), col=col, pch=points)
  p = anova(lm(rich~urb))$`Pr(>F)`[1]
  mtext(paste("p=",format(p, digits=2),sep=""), side=3, line=0)
  mtext(taxa, side = 3, line = 2, cex=1.5)
  pVal$pRichness[row] = p
  pVal$pRichnessWilcox[row] = wilcox.test(rich~urb, exact=F)$p.value
  
  ##inverse simpson
  inv = diversity(table[,-(1:3)], index="invsimpson")
  boxplot(inv~urb, ylab="Inverse Simpson Index")
  points(inv~jitter(x), col=col, pch=points)
  p = anova(lm(inv~urb))$`Pr(>F)`[1]
  mtext(paste("p=",format(p, digits=2),sep=""), side=3, line=0)
  pVal$pInvSimpson[row] = p
  pVal$pInvSimpsonWilcox[row] = wilcox.test(inv~urb, exact=F)$p.value
  
  ##shannon
  shan = diversity(table[,-(1:3)], index="shannon")
  boxplot(shan~urb, ylab="Shannon Index")
  points(shan~jitter(x), col=col, pch=points)
  p = anova(lm(shan~urb))$`Pr(>F)`[1]
  mtext(paste("p=",format(p, digits=2),sep=""), side=3, line=0)
  pVal$pShannon[row] = p
  pVal$pShannonWilcox[row] = wilcox.test(shan~urb, exact=F)$p.value
  
  ##evenness
  even = shan/log(rich)
  boxplot(even~urb, ylab="Evenness")
  points(even~jitter(x), col=col, pch=points)
  p = anova(lm(even~urb))$`Pr(>F)`[1]
  mtext(paste("p=",format(p, digits=2),sep=""), side=3, line=0)
  pVal$pShanEvenness[row] = p
  pVal$pShanEvennessWilcox[row] = wilcox.test(even~urb, exact=F)$p.value
  
  results = data.frame(richness=rich, invSimpson = inv, shannon=shan, evenness=even)
  results = cbind(table[,1:2], results)
  write.table(results, outputName,
              col.names = T, row.names = F, quote=F, sep="\t")
  
  return(pVal)
}

####logged data
###plot
pVal = data.frame(level=taxaLevels, 
                  pRichness = rep(NA, length(taxaLevels)), pRichnessWilcox = rep(NA, length(taxaLevels)), 
                  pInvSimpson = rep(NA, length(taxaLevels)), pInvSimpsonWilcox = rep(NA, length(taxaLevels)),
                  pShannon = rep(NA, length(taxaLevels)), pShannonWilcox = rep(NA, length(taxaLevels)),
                  pShanEvenness = rep(NA, length(taxaLevels)), pShanEvennessWilcox = rep(NA, length(taxaLevels)),
                  stringsAsFactors = F)
jpeg("minikraken_diversity_log.jpg", height=4200, width=4375, res=300)
layout(matrix(1:28, nrow=4, ncol=7, byrow=F), heights=1, widths=1)
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.9, oma=c(2,0,0,0), mar=c(4.1, 6, 4,.1))
for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  print(taxa)
  file = paste("minikraken_merged_taxaAsCol_logNorm_", taxa, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  pVal = drawPlots(table, taxa, paste("minikraken_diversity_log_", taxa, ".txt", sep=""), pVal)
}
##legend
par(oma=c(0,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", cex=1.5, horiz=T,
       legend=c("rural", "urban"),
       col=c("blue", "red"),
       pch=16)
dev.off()

###p-values
pVal$pAdjRichness = p.adjust(pVal$pRichness, method="BH")
pVal$pAdjRichnessWilcox = p.adjust(pVal$pRichnessWilcox, method="BH")
pVal$pAdjInvSimpson = p.adjust(pVal$pInvSimpson, method="BH")
pVal$pAdjInvSimpsonWilcox = p.adjust(pVal$pInvSimpsonWilcox, method="BH")
pVal$pAdjShannon = p.adjust(pVal$pShannon, method="BH")
pVal$pAdjShannonWilcox = p.adjust(pVal$pShannonWilcox, method="BH")
pVal$pAdjShanEvenness = p.adjust(pVal$pShanEvenness, method="BH")
pVal$pAdjShanEvennessWilcox = p.adjust(pVal$pShanEvennessWilcox, method="BH")
write.table(pVal, "minikraken_diversity_log_pValues.txt", col.names = T, row.names = F, quote=F, sep="\t")

####unlogged data
###plot
pVal = data.frame(level=taxaLevels, 
                  pRichness = rep(NA, length(taxaLevels)), pRichnessWilcox = rep(NA, length(taxaLevels)), 
                  pInvSimpson = rep(NA, length(taxaLevels)), pInvSimpsonWilcox = rep(NA, length(taxaLevels)),
                  pShannon = rep(NA, length(taxaLevels)), pShannonWilcox = rep(NA, length(taxaLevels)),
                  pShanEvenness = rep(NA, length(taxaLevels)), pShanEvennessWilcox = rep(NA, length(taxaLevels)),
                  stringsAsFactors = F)
jpeg("minikraken_diversity_unlog.jpg", height=4200, width=4375, res=300)
layout(matrix(1:28, nrow=4, ncol=7, byrow=F), heights=1, widths=1)
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.9, oma=c(2,0,0,0), mar=c(4.1, 6, 4,.1))
for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  print(taxa)
  file = paste("minikraken_merged_taxaAsCol_", taxa, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  pVal = drawPlots(table, taxa, paste("minikraken_diversity_unlog_", taxa, ".txt", sep=""), pVal)
}
##legend
par(oma=c(0,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", cex=1.5, horiz=T,
       legend=c("rural", "urban"),
       col=c("blue", "red"),
       pch=16)
dev.off()

###p-values
pVal$pAdjRichness = p.adjust(pVal$pRichness, method="BH")
pVal$pAdjRichnessWilcox = p.adjust(pVal$pRichnessWilcox, method="BH")
pVal$pAdjInvSimpson = p.adjust(pVal$pInvSimpson, method="BH")
pVal$pAdjInvSimpsonWilcox = p.adjust(pVal$pInvSimpsonWilcox, method="BH")
pVal$pAdjShannon = p.adjust(pVal$pShannon, method="BH")
pVal$pAdjShannonWilcox = p.adjust(pVal$pShannonWilcox, method="BH")
pVal$pAdjShanEvenness = p.adjust(pVal$pShanEvenness, method="BH")
pVal$pAdjShanEvennessWilcox = p.adjust(pVal$pShanEvennessWilcox, method="BH")
write.table(pVal, "minikraken_diversity_unlog_pValues.txt", col.names = T, row.names = F, quote=F, sep="\t")
