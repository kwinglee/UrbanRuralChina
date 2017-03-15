##minikraken figure of PCoAs, diversity, p-values and domain boxplots

rm(list=ls())
setwd("")

library(vegan)
library(ape)
library(ROCR)

tiff("minikrakenFig.tif", res=350, height=3500, width=11666)

##layout will be:
##row1 (10 total): domain pcoa, domain diversity, phylum pcoa, phy diversity, class pcoa, class diversity, order pcoa, order diversity, fam pcoa, fam div
##row2 (10 total): genus pcoa, gen div, otu pcoa, otu div, p plot, p legend, archaea, bacteria, viruses, legend
layout(matrix(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10,11,11,12,13,13,14,15,15,16,16,17,17,18,18,19), 
              nrow=2, ncol=15, byrow=2), 
       heights=1, widths=1)
# layout.show(19)


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
  pcoa <- capscale(table[,-(1:3)]~1,distance="bray")
  
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
       col=colors, pch=shapes, cex=4, cex.lab=3)
  mtext(taxa, side=3, line=1, cex=2.5, adj=1)
  mtext(lab, side=3, line=1.4, cex=2.5, adj=0)
}

##function to draw boxplots of Shannon diversity from given table
drawBoxplot <- function(table, taxa) {
  ##get diversity measures
  shan = diversity(table[,-(1:3)], index="shannon")
  
  ##plot urban rural
  par(mar=c(4.2, 4.8, 4, .1))
  boxplot(shan~table$ruralUrban, ylab="Shannon Index", xaxt="n", cex.lab=2.5, outline=F, ylim=range(shan))
  ##rotate x-axis to fit
  axis(1, at=1:2, labels=F)
  text(x=1:2, y=par("usr")[3]-.03, srt=25, adj=1, xpd=T, labels=c("rural", "urban"), cex=3)
  
  ##points
  colors = getColors(table)
  shapes = getShapes(table)
  ##jitter to space out
  x = ifelse(table$ruralUrban=="rural", 1, 2)
  points(x=jitter(x), y=shan, col=colors, pch=shapes, cex=2.2)
}

taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")
lab = c("A", "B", "C", "D", "E", "F", "G")

for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  file = paste("minikraken_merged_taxaAsCol_logNorm_", taxa, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  drawPcoa(table, taxa, lab[i])
  
  file = paste("minikraken_merged_taxaAsCol_", taxa, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  drawBoxplot(table, taxa)
}

############
###p value plot
xvals = 1:length(taxaLevels)
colors = c("black", "black")
lty = c(1, 4)
par(mar=plot.mar)

ymax = 3.2
ymin = 1

##polygon of significance region
plot.new()
plot.window(xlim=c(1,length(xvals)), ylim=c(ymin,ymax))
polygon(x=c(0, 0, 8, 8), y=c(-log10(0.05), ymax+.15, ymax+.15, -log10(0.05)), col="ivory2", border=NA)
par(new=T)
plot(1, type="n", xaxt="n", xlim=c(1,length(xvals)), ylim=c(0,ymax), ylab="-log10(urban/rural P-value)", xlab="", cex.lab=2.4)
axis(1, at=1:length(xvals), labels=F) #add x-axis ticks
text(x=1:length(xvals), y=par("usr")[3]-.1, srt=22, adj=1, xpd=T, labels=c("domain", "phylum", "class", "order", "family", "genus", "species"), cex=2.3) #rotate x-axis labels

mtext("H", side=3, line=1.4, cex=2, adj=0)


##function to draw lines of the given color for the given rural urban pvalues (pvals)
drawRuralUrbanLines <- function(pvals, color, lty) {
  lines(x=xvals, y=-log10(pvals), lty=lty, col=color, lwd=5)
  points(x=xvals, y=-log10(pvals), pch=15, col=color, cex=2.5)
}

####MDS p values
mds1urbanRural = rep(NA, length(taxaLevels))
for(i in xvals) {
  taxa = taxaLevels[i]
  table = read.table(paste("minikraken_pcoaModel_pValues_", taxa, ".txt", sep=""), 
                     sep="\t", header=T, colClasses=c("character", rep("numeric", 8)))
  mds1row = table$names=="MDS1"
  mds1urbanRural[i] = table$pValuesUrbanRural[mds1row]
}
##MDS1 lines
drawRuralUrbanLines(mds1urbanRural, colors[1], lty[1])

####Shannon
table = read.table(paste("minikraken_diversity_unlog_pValues.txt", sep=""),  
                   header=T, sep="\t", colClasses = c("character", rep("numeric", 16)))
drawRuralUrbanLines(table$pShannon, colors[2], lty[2])

#######################
####boxplots of domains
table = read.table("minikraken_merged_taxaAsCol_logNorm_domain.txt", 
                   header=T, sep="\t", 
                   colClasses=c("character", "numeric", "character", rep("numeric", 3)))
colors = getColors(table)
shapes = getShapes(table)
rurUrb = factor(table$ruralUrban)
x = ifelse(table$ruralUrban=="rural", 1, 2)
pValues = read.table("minikraken_otuModel_pValues_domain.txt",
                     header=T, sep="\t", colClasses = c("character", rep("numeric", 9)))

##for the given p-value, return the correct number of stars
getStars <- function(pval) {
  if(pval < 0.001) {
    return("***")
  } else if(pval >= 0.001 && pval < 0.01) {
    return("**")
  } else if(pval >= 0.01 && pval < 0.05) {
    return("*")
  } else if(pval >= 0.05 && pval < 0.1) {
    return(".")
  } else {
    return("")
  }
}

##for the given relative abundance, draw a boxplot with the given title and pvalue
drawRelAbunBoxplot <- function(relAbun, title, pval) {
  par(mar=c(4.2, 7, 4, 3))
  boxplot(relAbun~rurUrb, ylab="log relative abundance", xaxt="n", 
          cex.lab=2.5, outline=F, ylim=range(relAbun))
  mtext(title, side=3, line=1, cex=2.5)
  ##rotate x-axis to fit
  axis(1, at=1:2, labels=F, cex.axis=3, line=0)
  axis(1, at=1:2, labels=c("rural", "urban"), cex.axis=3, line=.5, tick=F)
  ##points
  ##jitter to space out
  points(x=jitter(x), y=relAbun, col=colors, pch=shapes, cex=3)
  ##p value
  mtext(getStars(pval), side=1, line=1, cex=5)
}

drawRelAbunBoxplot(table$Archaea, "Archaea", pValues$adjustedPurbanRural[pValues$names=="Archaea"])
mtext("I", side=3, line=1.4, cex=2, adj=0)
drawRelAbunBoxplot(table$Bacteria, "Bacteria", pValues$adjustedPurbanRural[pValues$names=="Bacteria"])
drawRelAbunBoxplot(table$Viruses, "Viruses", pValues$adjustedPurbanRural[pValues$names=="Viruses"])

##legend
par(mar=leg.mar)
plot(1, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
legend("top", 
       title = expression(underline("Key for A-G & I")),
       legend=c("rural", "urban"),
       col=c("blue", "red"),
       pch=16,
       cex=2.7)
legend("top", inset = c(0, .45), 
       title = "",
       legend=c("MDS1", "Shannon", "p<0.05"),
       col=c("black", "black", "ivory2"),
       lty=c(lty, 1), lwd=c(5, 5, 20), cex=2.4)
legend("top", inset = c(0, .45), 
       legend="",
       title = expression(underline("Key for H")),
       bty="n", cex=2.7)
dev.off()
