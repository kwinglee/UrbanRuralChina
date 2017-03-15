##Figure 1
##16S microbiome figure of pcoas, diversity, pvalues and roc

rm(list=ls())
setwd("")

library(vegan)
library(ape)
library(ROCR)

tiff("microbiomeFig.tif", res=350, height=3500, width=9333)

##layout will be:
##row1 (9 total): phylum pcoa, phy diversity, class pcoa, class diversity, order pcoa, order diversity, fam pcoa, fam div, legend
##row2 (8 total): genus pcoa, gen div, otu pcoa, otu div, p plot, p legend, roc, roc legend
layout(matrix(c(1:16, 16, 16), nrow=2, ncol=9, byrow=2), heights=1, 
       widths=c(rep(c(2, 1), 4), 1))

plot.mar = c(4.2, 7, 4,.1) #plot margins
leg.mar = c(4.2, .1, 4,.1) #legend margins
par(cex.axis=2, cex.lab=2)

##############
###pcoa and diversity plots
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
drawPcoa <- function(table, taxa, lab) {
  colors = getColors(table)
  shapes = getShapes(table)
  
  ##using vegan capscale with bray curtis distance
  pcoa <- capscale(table[,-(1:5)]~1,distance="bray")
  
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
       col=colors, pch=shapes, cex=3.3, cex.lab=2.5)
  mtext(taxa, side=3, line=1, cex=2.5, adj=1)
  mtext(lab, side=3, line=1.4, cex=2.5, adj=0)
}

##function to draw boxplots of Shannon diversity from given diversity table
drawBoxplot <- function(table) {
  ##get diversity measures
  shan = table$shannon
  
  ##plot urban rural
  par(mar=c(4.2, 4.8, 4, .1))
  boxplot(shan~table$ruralUrban, ylab="Shannon Index", xaxt="n", cex.lab=2.5)
  ##rotate x-axis to fit
  axis(1, at=1:2, labels=F)
  text(x=1:2, y=par("usr")[3]-.03, srt=25, adj=1, xpd=T, labels=c("rural", "urban"), cex=2.5)
  
  ##points
  colors = getColors(table)
  shapes = getShapes(table)
  ##jitter to space out
  x = ifelse(table$ruralUrban=="rural", 1, 2)
  points(x=jitter(x), y=shan, col=colors, pch=shapes, cex=2.5)
}

taxaLevels <- c("phylum","class","order","family","genus")
lab = c("A", "B", "C", "D", "E")

for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  fileName <- paste(taxa, "_taxaAsColumnsLogNorm_WithMetadata.txt", sep ="")
  table <-read.table(fileName,header=TRUE,sep="\t")
  numCols <- ncol(table)
  table <-read.table(fileName,header=TRUE,sep="\t",colClasses=c("character", "numeric", "numeric", "character", "character", rep("numeric", numCols-5)))
  table = table[table$readNumber==1,]
  
  drawPcoa(table, taxa, lab[i])
  
  div = read.table(paste("rarefyDiversity_diversityValues_", taxa, "_average.txt", sep=""),
                   header=T, sep="\t", colClasses = c("character", rep("numeric", 4)))
  div = merge(table[,1:5], div, by="sampleID")
  
  drawBoxplot(div)
  
  if(taxa == "family") { #time to draw legend
    par(mar=leg.mar)
    plot(1, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
    legend("topleft", c("rural", "urban", "timepoint 1", "timepoint 2"), 
           col=c("blue", "red", "gray", "gray"), pch=c(15, 15, 16, 17), cex=2.55)
  }
}

##OTU
table = read.table("abundantOTUForwardTaxaAsColumnsLogNormalWithMetadata.txt", sep="\t", header=T, colClasses=c("character", "numeric", "character", "character", rep("numeric", 870)))
table = cbind(table$sampleID, readNumber=rep(1, nrow(table)), table[,c(2:4, 8:ncol(table))]) #add read number, remove diversity
names(table)[1] = "sampleID"
table = table[,-c(6:8)] #remove diversity
drawPcoa(table, "OTU", "F")
div = read.table("rarefyDiversity_diversityValues_otu_average.txt",
                 header=T, sep="\t", colClasses = c("character", rep("numeric", 4)))
table$sampleID = paste(table$sampleID, "_1", sep="")
div = merge(table[,1:5], div, by="sampleID")
drawBoxplot(div)

############
###p value plot
xvals = 1:6
taxaLevels = c("phylum", "class", "order", "family", "genus", "otu")
colors = c("black", "black")
lty = c(1, 4)
par(mar=plot.mar)

ymax = 4.1

##polygon of significance region
plot.new()
plot.window(xlim=c(1,length(xvals)), ylim=c(0,ymax))
polygon(x=c(0, 0, 7, 7), y=c(-log10(0.05), ymax+.15, ymax+.15, -log10(0.05)), col="ivory2", border=NA)
par(new=T)
plot(1, type="n", xaxt="n", xlim=c(1,length(xvals)), ylim=c(0,ymax), ylab="-log10(urban/rural P-value)", xlab="", cex.lab=2.4)
axis(1, at=1:length(xvals), labels=F) #add x-axis ticks
text(x=1:length(xvals), y=par("usr")[3]-.1, srt=22, adj=1, xpd=T, labels=c("phylum", "class", "order", "family", "genus", "OTU"), cex=2.3) #rotate x-axis labels

mtext("G", side=3, line=1.4, cex=2, adj=0)


##function to draw lines of the given color for the given rural urban pvalues (pvals)
drawRuralUrbanLines <- function(pvals, color, lty) {
  lines(x=xvals, y=-log10(pvals), lty=lty, col=color, lwd=5)
  points(x=xvals, y=-log10(pvals), pch=15, col=color, cex=2.5)
}

####MDS p values
mds1urbanRural = rep(NA, length(taxaLevels))
for(i in xvals) {
  taxa = taxaLevels[i]
  table = read.table(paste("pcoaModel_pValues_", taxa, ".txt", sep=""), sep="\t", header=T, colClasses=c("character", rep("numeric", 12)))
  mds1row = table$names=="MDS1"
  mds1urbanRural[i] = table$pValuesUrbanRural[mds1row]
}
##MDS1 lines
drawRuralUrbanLines(mds1urbanRural, colors[1], lty[1])

##legend
par(mar=leg.mar)
plot(1, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
legend("topleft",
       legend=c("MDS1", "p<0.05"),
       col=c("black", "ivory2"),
       lty=c(lty[1], 1), lwd=c(5, 20), cex=2.55)

#######################
####ROC
taxaLevels = c("phylum", "class", "order", "family", "genus", "otu")
colors = c("blue", "black", "red", "purple", "turquoise", "darkgreen")

##set up plot
par(mar=c(plot.mar[1:3], .2))
plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate", ylab="True positive rate", cex.lab=2.4)
mtext("H", side=3, line=1.5, cex=2, adj=0)

##predictions
auc = rep(NA, length(taxaLevels))
for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  table = read.table(paste("cforest_", taxa, ".txt", sep=""), sep="\t", header=T, colClasses="numeric")
  pred = prediction(table$predicted, table$actual)
  perf = performance(pred, measure="tpr", x.measure="fpr")
  lines(x=unlist(slot(perf, "x.values")), y=unlist(slot(perf, "y.values")), col=colors[i], lwd=5)
  a = performance(pred,"auc")
  auc[i] = unlist(slot(a, "y.values"))
}

##legend
legend("bottomright",
       legend = c("Input Data", "OTU", "Genus", "Family", "Order", "Class", "Phylum                "),
       col=c("white", rev(colors)), lty=1, lwd=2, cex=2.35)
text(x=.97, y=seq(0.56, 0.04, length.out=length(auc)+1), labels=c("AUC", rev(format(auc, digits=2))), cex=2.35)

dev.off()
