##plot the minikraken diversity measures

setwd("")

library(vegan)
library(nlme)

taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")

##get metadata
meta = read.table("phylum_taxaAsColumnsLogNorm_withMetadata.txt", header=T, sep="\t", colClasses = "character")
meta = meta[,1:5]
meta = meta[meta$readNumber=="1",]
meta$sampleID = gsub("_1", "", meta$sampleID)

##p-values
pVal = read.table("minikraken_diversity_unlog_pValues.txt", header=T, sep="\t",
                  colClasses = c("character", rep("numeric", 16)))

tiff("minikrakenDivSuppFig_unlog.tif", height=3200, width=4260, res=300)
layout(matrix(1:21, nrow=3, ncol=7, byrow=F), heights=1, widths=1)
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.9, oma=c(2,0,0,.3), mar=c(4.1, 6, 4,.1))

for(lev in taxaLevels) {
  print(lev)
  table = read.table(paste("minikraken_diversity_unlog_", lev, ".txt", sep=""), header=T, sep="\t",
                     colClasses = c("character", rep("numeric", 5)))
  table = merge(meta, table, by="sampleID")
  
  ##plot
  urb = factor(table$ruralUrban)
  x=ifelse(urb=="rural", 1, 2)
  col=ifelse(urb=="rural", "blue", "red")
  points=ifelse(table$timepoint=="first_A", 16, 17)
  
  boxplot(table$invSimpson~urb, ylab="Inverse Simpson Index", outline=F, ylim=range(table$invSimpson))
  points(table$invSimpson~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pVal$pInvSimpson[pVal$level==lev], digits=2),sep=""), side=3, line=0)
  if(lev == "otu") {
    lev = "OTU"
  }
  mtext(lev, side = 3, line = 2, cex=1.5)
  if(lev == "domain") {
    mtext("A", side = 2, line = 4, cex=2, padj=0, las=1, at=1.031)
  }
  
  boxplot(table$richness~urb, ylab="Richness", outline=F, ylim=range(table$richness)) 
  points(table$richness~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pVal$pRichness[pVal$level==lev], digits=2),sep=""), side=3, line=0)
  if(lev == "domain") {
    mtext("B", side = 2, line = 4, cex=2, padj=0, las=1, at=4.6)
  }
  
  boxplot(table$evenness~urb, ylab="Evenness", outline=F, ylim=range(table$evenness))
  points(table$evenness~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pVal$pShanEvenness[pVal$level==lev], digits=2),sep=""), side=3, line=0)
  if(lev == "domain") {
    mtext("C", side = 2, line = 4, cex=2, padj=0, las=1, at=.075)
  }
}

##legend
par(oma=c(.2,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", cex=1.5, horiz=T,
       legend=c("rural", "urban"),
       col=c("blue", "red"),
       pch=c(16))
dev.off()