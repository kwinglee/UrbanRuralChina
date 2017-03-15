##make separate figures for each level of pcoa + logged diversity
##output tables of PCoA site scores

rm(list=ls())
library("vegan")
library("ape")

setwd("")

taxaLevels <- c("phylum","class","order","family","genus")

##function to get get colors
getColors <- function(table) {
  colors = rep("blue", nrow(table)) #rural
  colors[table$ruralUrban=="urban"] = "red"
  return(colors)
}

##function to get point shapes
getShapes <- function(table) {
  shapes = rep(16, nrow(table)) #first time point
  shapes[grepl("B", table$sampleID, fixed=F)] = 17
  return(shapes)
}

##function to run the PCOA on the given table, writing results using the taxa name
##writes the corrected pcoa aces
drawPcoa <- function(table, taxa) {
  colors = getColors(table)
  shapes = getShapes(table)
  
  ##using vegan capscale
  pcoa <- capscale(table[,-(1:5)]~1,distance="bray")
  plot(x=pcoa$CA$u[,1], y=pcoa$CA$u[,2], xlab="MDS1", ylab="MDS2", main=taxa, col=colors, pch=shapes)

  ##write axes
  combine = cbind(table[,1:5], pcoa$CA$u)
  write.table(combine, sep="\t", file=paste("pcoaCorrected_", taxa, ".txt",sep=""), quote=F, row.names=F, col.names=T)
  
}

##function to draw boxplots of Shannon diversity from given table
##returns the p-value
drawBoxplot <- function(table, taxa) {
  ##get diversity measures
  shan = diversity(table[,-(1:5)], index="shannon")
  
  ##t tests
  p.shan = t.test(shan~table$ruralUrban)$p.value
  p.time = t.test(shan~table$timepoint)$p.value
  
  ##plot urban rural
  boxplot(shan~table$ruralUrban, main=paste("Urban vs. Rural\nP=", round(p.shan, digits=3), sep=""), ylab="Shannon Index", main.cex=0.5)
  
  ##points
  colors = getColors(table)
  shapes = getShapes(table)
  points(x=factor(table$ruralUrban), y=shan, col=colors, pch=shapes)
  
  ##plot time
  boxplot(shan~table$timepoint, main=paste("Timepoint\nP=", round(p.time, digits=3), sep=""), ylab="Shannon Index", main.cex=0.5, xaxt="n")
  axis(1, at=1:2, labels=c("1", "2"))
  points(x=factor(table$timepoint), y=shan, col=colors, pch=shapes)
  
  return(c(p.shan, p.time))
}

h=1200
w=1800

pvals = data.frame(level=character(), pShannonUrbanRural=numeric(), pShannonTime=numeric())

for(taxa in taxaLevels ) {
  fileName <- paste(taxa, "_taxaAsColumnsLogNorm_WithMetadata.txt", sep ="")
  table <-read.table(fileName,header=TRUE,sep="\t")
  numCols <- ncol(table)
  table <-read.table(fileName,header=TRUE,sep="\t",colClasses=c("character", "numeric", "numeric", "character", "character", rep("numeric", numCols-5)))
  table = table[table$readNumber==1,]
  
  jpeg(paste("pcoaWithDiversity_", taxa, ".jpg", sep=""), res=200, height=h, width=w)
  # layout(matrix(c(1,2), nrow=1, ncol=2, byrow=T), widths=c(2,1), heights=c(1,1))
  layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow=T), widths=c(2,1,1), heights=c(1,1,1))
  drawPcoa(table, taxa)
  p = drawBoxplot(table, taxa)
  dev.off()
  
  pvals = rbind(pvals, data.frame(level=taxa, pShannonUrbanRural=p[1], pShannonTime=p[2]))
}

##OTU
table = read.table("abundantOTUForwardTaxaAsColumnsLogNormalWithMetadata.txt", sep="\t", header=T, colClasses=c("character", "numeric", "character", "character", rep("numeric", 870)))
table = cbind(table$sampleID, readNumber=rep(1, nrow(table)), table[,2:ncol(table)]) #add read number
names(table)[1] = "sampleID"
table = table[,-c(6:8)] #remove diversity
jpeg("pcoaWithDiversity_otu.jpg", res=200, height=h, width=w)
layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow=T), widths=c(2,1,1), heights=c(1,1,1))
drawPcoa(table, "OTU")
p = drawBoxplot(table, "OTU")
dev.off()
pvals = rbind(pvals, data.frame(level="OTU", pShannonUrbanRural=p[1], pShannonTime=p[2]))
