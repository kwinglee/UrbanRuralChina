##look at diversity measures for HUMAnN using vegan (vs. using the HUMAnN output)

rm(list=ls())

library(vegan)

setwd("")

levels = c("module", "pathway")

##draws the plots for richness, inverse simpson, shannon and shannon evenness for the given table
##uses the given taxa name and writes the measures to the file outputName, puts the values in pVal
drawPlots <- function(table, taxa, outputName, pVal) {
  urb = factor(table$ruralUrban)
  col=ifelse(urb=="rural", "blue", "red")
  points=16
  x=ifelse(urb=="rural", 1, 2)
  row = pVal$level==taxa
  
  ##richness
  rich = specnumber(table[,-(1:6)])
  boxplot(rich~urb, ylab="Richness")
  points(rich~jitter(x), col=col, pch=points)
  p = anova(lm(rich~urb))$`Pr(>F)`[1]
  mtext(paste("p=",format(p, digits=2),sep=""), side=3, line=0)
  mtext(taxa, side = 3, line = 2, cex=1.5)
  pVal$pRichness[row] = p
  pVal$pRichnessWilcox[row] = wilcox.test(rich~urb, exact=F)$p.value
  
  ##inverse simpson
  inv = diversity(table[,-(1:6)], index="invsimpson")
  boxplot(inv~urb, ylab="Inverse Simpson Index")
  points(inv~jitter(x), col=col, pch=points)
  p = anova(lm(inv~urb))$`Pr(>F)`[1]
  mtext(paste("p=",format(p, digits=2),sep=""), side=3, line=0)
  pVal$pInvSimpson[row] = p
  pVal$pInvSimpsonWilcox[row] = wilcox.test(inv~urb, exact=F)$p.value
  
  ##shannon
  shan = diversity(table[,-(1:6)], index="shannon")
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
pVal = data.frame(level=levels, 
                  pRichness = rep(NA, length(levels)), pRichnessWilcox = rep(NA, length(levels)), 
                  pInvSimpson = rep(NA, length(levels)), pInvSimpsonWilcox = rep(NA, length(levels)),
                  pShannon = rep(NA, length(levels)), pShannonWilcox = rep(NA, length(levels)),
                  pShanEvenness = rep(NA, length(levels)), pShanEvennessWilcox = rep(NA, length(levels)),
                  stringsAsFactors = F)
jpeg("humann_vegan_diversity_log.jpg", height=4200, width=1200, res=300)
layout(matrix(1:8, nrow=4, ncol=2, byrow=F), heights=1, widths=1)
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.9, oma=c(2,0,0,0), mar=c(4.1, 6, 4,.1))
for(i in 1:length(levels)) {
  lev = levels[i]
  print(lev)
  file = paste("humann_keggAsCol_log_", lev, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", ncol-2)))
  
  pVal = drawPlots(table, lev, paste("humann_vegan_diversity_log_", lev, ".txt", sep=""), pVal)
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
write.table(pVal, "humann_vegan_diversity_log_pValues.txt", col.names = T, row.names = F, quote=F, sep="\t")

####unlogged data
###plot
pVal = data.frame(level=levels, 
                  pRichness = rep(NA, length(levels)), pRichnessWilcox = rep(NA, length(levels)), 
                  pInvSimpson = rep(NA, length(levels)), pInvSimpsonWilcox = rep(NA, length(levels)),
                  pShannon = rep(NA, length(levels)), pShannonWilcox = rep(NA, length(levels)),
                  pShanEvenness = rep(NA, length(levels)), pShanEvennessWilcox = rep(NA, length(levels)),
                  stringsAsFactors = F)
jpeg("humann_vegan_diversity_unlog.jpg", height=4200, width=1200, res=300)
layout(matrix(1:8, nrow=4, ncol=2, byrow=F), heights=1, widths=1)
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.9, oma=c(2,0,0,0), mar=c(4.1, 6, 4,.1))
for(i in 1:length(levels)) {
  lev = levels[i]
  print(lev)
  file = paste("humann_keggAsCol_withRurUrb_", lev, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", ncol-2)))
  
  pVal = drawPlots(table, lev, paste("humann_vegan_diversity_unlog_", lev, ".txt", sep=""), pVal)
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
write.table(pVal, "humann_vegan_diversity_unlog_pValues.txt", col.names = T, row.names = F, quote=F, sep="\t")
