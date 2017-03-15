##plot the rarefied 16S diversity measures

setwd("")

library(vegan)
library(nlme)

levels = c("phylum","class","order","family","genus", "otu")

##get metadata
meta = read.table("phylum_taxaAsColumnsLogNorm_withMetadata.txt", header=T, sep="\t", colClasses = "character")
meta = meta[,1:5]

##get p-values
pRich = read.table("rarefyDiversity_pValues_richness.txt", header=T, sep="\t", 
                   colClasses = c("character", rep("numeric", 6)))
pInvSimp = read.table("rarefyDiversity_pValues_invsimpson.txt", header=T, sep="\t", 
                   colClasses = c("character", rep("numeric", 6)))
pEven = read.table("rarefyDiversity_pValues_evenness.txt", header=T, sep="\t", 
                   colClasses = c("character", rep("numeric", 6)))

tiff("microbiome16sDivSuppFig_unlogRarefy.tif", height=3200, width=3650, res=300)
layout(matrix(1:18, nrow=3, ncol=6, byrow=F), heights=1, widths=1)
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.9, oma=c(2,0,0,.3), mar=c(4.1, 6, 4,.1))

for(lev in levels) {
  print(lev)
  table = read.table(paste("rarefyDiversity_diversityValues_", lev, "_average.txt", sep=""), header=T, sep="\t")
  table = merge(meta, table, by="sampleID")
  
  ##plot
  urb = factor(table$ruralUrban)
  x=ifelse(urb=="rural", 1, 2)
  col=ifelse(urb=="rural", "blue", "red")
  points=ifelse(table$timepoint=="first_A", 16, 17)
  
  boxplot(table$invSimpson~urb, ylab="Inverse Simpson Index", outline=F, ylim=range(table$invSimpson))
  points(table$invSimpson~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pInvSimp$pUrbanRural[pInvSimp$level==lev], digits=2),sep=""), side=3, line=0)
  if(lev == "otu") {
    mtext("OTU", side = 3, line = 2, cex=1.5)
  } else {
    mtext(lev, side = 3, line = 2, cex=1.5)
  }
  if(lev == "phylum") {
    mtext("A", side = 2, line = 4, cex=2, padj=0, las=1, at=3.2)
  }
  
  boxplot(table$richness~urb, ylab="Richness", outline=F, ylim=range(table$richness)) 
  points(table$richness~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pRich$pUrbanRural[pRich$level==lev], digits=2),sep=""), side=3, line=0)
  if(lev == "phylum") {
    mtext("B", side = 2, line = 4, cex=2, padj=0, las=1, at=11)
  }
  
  boxplot(table$evenness~urb, ylab="Evenness", outline=F, ylim=range(table$evenness))
  points(table$evenness~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pEven$pUrbanRural[pEven$level==lev], digits=2),sep=""), side=3, line=0)
  if(lev == "phylum") {
    mtext("C", side = 2, line = 4, cex=2, padj=0, las=1, at=.7)
  }
}

##legend
par(oma=c(.2,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", cex=1.5, horiz=T,
       legend=c("rural", "urban", "timepoint 1", "timepoint 2"),
       col=c("blue", "red", "gray", "gray"),
       pch=c(15, 15, 16, 17))
dev.off()