##plot the other HUMAnN diversity measures

setwd("")

levels = c("module", "pathway")

##get p-values
pVal = read.table("humann_vegan_diversity_unlog_pValues.txt", header=T, sep="\t", 
                   colClasses = c("character", rep("numeric", 16)))

tiff("humannDivSuppFig_unlog.tif", height=3200, width=1200, res=300)
par(cex.axis=1.5, cex.lab=1.5, cex.main=1.9, oma=c(3.1,0,0,.4), mar=c(4.1, 6, 4,.1),
    mfcol=c(3,2))

for(lev in levels) {
  print(lev)
  table = read.table(paste("humann_vegan_diversity_unlog_", lev, ".txt", sep=""), header=T, sep="\t", 
                     colClasses=c(rep("character", 2), rep("numeric", 4)))
  
  ##plot
  urb = factor(table$ruralUrban)
  x=ifelse(urb=="rural", 1, 2)
  col=ifelse(urb=="rural", "blue", "red")
  points=16
  
  boxplot(table$invSimpson~urb, ylab="Inverse Simpson Index", outline=F, ylim=range(table$invSimpson))
  points(table$invSimpson~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pVal$pInvSimpson[pVal$level==lev], digits=2),sep=""), side=3, line=0)
  mtext(lev, side = 3, line = 2, cex=1.5)
  if(lev == "module") {
    mtext("A", side = 2, line = 4, cex=2, padj=0, las=1, at=72.5)
  }
  
  boxplot(table$richness~urb, ylab="Richness", outline=F, ylim=range(table$richness)) 
  points(table$richness~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pVal$pRichness[pVal$level==lev], digits=2),sep=""), side=3, line=0)
  if(lev == "module") {
    mtext("B", side = 2, line = 4, cex=2, padj=0, las=1, at=208)
  }
  
  boxplot(table$evenness~urb, ylab="Evenness", outline=F, ylim=range(table$evenness))
  points(table$evenness~jitter(x), col=col, pch=points)
  mtext(paste("p=",format(pVal$pShanEvenness[pVal$level==lev], digits=2),sep=""), side=3, line=0)
  if(lev == "module") {
    mtext("C", side = 2, line = 4, cex=2, padj=0, las=1, at=.835)
  }
}

##legend
par(oma=c(.01,0,0,0), mar=c(0.3,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", cex=1.5, ncol=2,
       legend=c("rural", "urban", "timepoint 1", "timepoint 2"),
       col=c("blue", "red", "gray", "gray"),
       pch=c(15, 15, 16, 17))
dev.off()