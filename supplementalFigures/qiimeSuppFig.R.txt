##draw QIIME PCoAs for supplemental figure

rm(list=ls())
setwd("")

library("vegan")
library(ROCR)

tiff("qiimeSuppFig.tif", res=200, height=2000, width=3500)
##layout:
##row1: phylum, class, order, legend
##row2: family, genus, otu, blank
layout(matrix(c(1:8), nrow=2, ncol=4, byrow=2), heights=1, widths=rep(c(rep(1, 3), .2), 2))
plot.mar = c(4, 7, 4,.1) #plot margins
leg.mar = c(4, .1, 4,.1) #legend margins
par(cex.axis=1.5, cex.lab=1.5, cex.main=2.5)

##function to get get colors
getColors <- function(table) {
  colors = rep("blue", nrow(table)) #rural
  colors[table$ruralUrban=="urban"] = "red"
  return(colors)
}

##function to run the PCOA on the given table, writing results using the taxa name
drawPcoa <- function(table, taxa, lab) {
  colors = getColors(table)
  
  ##using vegan capscale
  pcoa <- capscale(table[,-(1:2)]~1,distance="bray")
  ##get percent variance
  eig = eigenvals(pcoa)
  if(any(eig<0)) {#need to check there are no negatives
    warning(paste("NEGATIVE EIGENVALUES-percent variance is incorrect for ", taxa, sep=""))
  }
  var = eig/sum(eig)*100
  par(mar=plot.mar)
  plot(x=pcoa$CA$u[,1], y=pcoa$CA$u[,2], 
       xlab=paste("MDS1 (", format(var[[1]], digits=0), "%)", sep=""), 
       ylab=paste("MDS2 (", format(var[[2]], digits=0), "%)", sep=""), 
       pch=16, xlim=c(-.45, .45), ylim=c(-.4, .4), main=taxa, cex=2, col=colors)
  mtext(lab, side=3, line=1.5, cex=2, adj=0)
}

taxaLevels <- c("phylum","class","order","family","genus", "OTU")
lab = c("A", "B", "C", "D", "E", "F")
fileExts = c("_L2", "_L3", "_L4", "_L5", "_L6", "")

for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  print(taxa)
  fileName <- paste("qiime_otu_table", fileExts[i], "_taxaAsCol_logNorm_with_metadata.txt", sep ="")
  table <-read.table(fileName,header=TRUE,sep="\t")
  numCols <- ncol(table)
  table <-read.table(fileName,header=TRUE,sep="\t",colClasses=c("character", "character", rep("numeric", numCols-2)))

  drawPcoa(table, taxa, lab[i])
  
  if(taxa == "order") { #time to draw legend
    par(mar=leg.mar)
    plot(1, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
    legend("topleft", c("rural", "urban"), 
           col=c("blue", "red"), pch=16, cex=2)
  }
}

dev.off()
