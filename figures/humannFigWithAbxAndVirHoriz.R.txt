##humann figure of unlogged PCoAs, diversity, and pvalues, with antibiotic and virulence plots

rm(list=ls())
setwd(")

library(vegan)
library(ape)
library(ROCR)

tiff("humannFigWithAbxAndVirHoriz.tif", res=350, height=3500, width=4667)

##layout will be:
##row1 (7 total): module pcoa, module div, pathway pcoa, pathway div, p plot, antibiotic plot, legend
##alt row1 (3 total): module pcoa, module div, p plot, legend
##alt row2 (4 total): pathway pcoa, pathway div, abx plot, vir plot
layout(matrix(c(1,1,2,3,4,5,5,6,7,8), nrow=2, ncol=5, byrow=2), heights=1, 
       widths=1)
# layout.show(8)

plot.mar = c(4.2, 7, 4,.1) #plot margins
leg.mar = c(4.2, 2, 4,.1) #legend margins
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

##function to run the PCOA on the given table, with the level name and ABC label lab
drawPcoa <- function(table, level, lab) {
  colors = getColors(table)
  shapes = getShapes(table)
  
  ##using vegan capscale with bray curtis distance
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
       col=colors, pch=shapes, cex=4, cex.lab=2.5)
  mtext(level, side=3, line=1, cex=2.5, adj=1)
  mtext(lab, side=3, line=1.4, cex=2.5, adj=0)
}

##function to draw boxplots of Shannon diversity from given table
drawBoxplot <- function(table) {
  ##get diversity measures
  shan = diversity(table[,-(1:2)], index="shannon")
  
  ##plot urban rural
  par(mar=c(4.2, 5, 4, .1))
  boxplot(shan~table$ruralUrban, ylab="Shannon Index",  cex.lab=2.9, 
          ylim=range(shan), outline = F, xaxt="n")
  axis(side=1, at=c(1,2), labels=c("rural", "urban"), line=0, cex.axis=2.3)
 
  ##points
  colors = getColors(table)
  shapes = getShapes(table)
  ##jitter to space out
  x = ifelse(table$ruralUrban=="rural", 1, 2)
  points(x=jitter(x), y=shan, col=colors, pch=shapes, cex=2.7)
}

#####module
table = read.table("humann_keggAsCol_withRurUrb_module.txt", header=T, sep="\t")
table = table[,-(3:6)]
drawPcoa(table, "module", "A")
drawBoxplot(table)

############
###p value plot
levels = c("module", "pathway")
xvals = 1:length(levels)
colors = c("black", "black")
lty = c(1, 4)
par(mar=plot.mar)

ymax = 4.1

##polygon of significance region
plot.new()
plot.window(xlim=c(.5,2.5), ylim=c(1,ymax))
polygon(x=c(0, 0, 7, 7), y=c(-log10(0.05), ymax+.15, ymax+.15, -log10(0.05)), col="ivory2", border=NA)
par(new=T)
plot(1, type="n", xaxt="n", xlim=c(.9, 2.1), ylim=c(1,ymax), ylab="-log10(urban/rural P-value)", xlab="", cex.lab=2.4)
axis(1, at=1:length(xvals), labels=F) #add x-axis ticks
axis(side=1, at=c(1,2), labels=c("module", "pathway"), line=0, cex.axis=2.4)

mtext("C", side = 2, line = 2.2, cex=2.5, padj=0, las=1, at=4.41)


##function to draw lines of the given color for the given rural urban pvalues (pvals)
drawRuralUrbanLines <- function(pvals, color, lty) {
  lines(x=xvals, y=-log10(pvals), lty=lty, col=color, lwd=6)
  points(x=xvals, y=-log10(pvals), pch=15, col=color, cex=2.5)
}

####MDS p values
mds1urbanRural = rep(NA, length(levels))
for(i in xvals) {
  lev = levels[i]
  table = read.table(paste("humann_pcoaModel_pValues_unlog_", lev, ".txt", sep=""), sep="\t", header=T, 
                     colClasses=c("character", rep("numeric", 8)))
  mds1row = table$names=="MDS1"
  mds1urbanRural[i] = table$pValuesUrbanRural[mds1row]
}
##MDS1 lines
# drawRuralUrbanLines(p.adjust(mds1urbanRural, "BH"), colors[1], lty[1])
drawRuralUrbanLines(mds1urbanRural, colors[1], lty[1])

####Shannon
table = read.table("humann_vegan_diversity_unlog_pValues.txt", header=T, sep="\t", 
                   colClasses = c("character", rep("numeric", 16)))
drawRuralUrbanLines(table$pShannon, colors[2], lty[2])

#########
####legend
par(mar=c(.1, .1, .3, .1))
leg.cex = 2.9
plot(1, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
legend(x=.5, y=.7, xjust=.5, yjust=.5, bty="n",
       title = expression(underline("Key for A, B, D, E")),
       legend=c("rural", "urban"),
       col=c("blue", "red"),
       pch=16, cex=leg.cex)
legend(x=.5, y=.5, xjust=.5, bty="n",
       title = expression(underline("Key for C")),
       legend=c("MDS1", "Shannon", "p<0.05"),
       col=c("black", "black", "ivory2"),
       lty=c(lty, 1), lwd=c(5, 5, 20), cex=leg.cex)


##########
#####pathway
table = read.table("humann_keggAsCol_withRurUrb_pathway.txt", header=T, sep="\t")
table = table[,-(3:6)]
drawPcoa(table, "pathway", "B")
drawBoxplot(table)

############
####antibiotic plot
table = read.table("pro_homolog_results.txt", header=T, sep="\t", quote="",
                   colClasses=c("character", rep("numeric", 1411)))
names(table)[1] = "sampleID"
##add metadata
meta = read.table("genus_taxaAsColumns_relAbunFwd.txt", sep="\t", header=T,
                  colClasses=c(rep("character", 3), rep("numeric", 347)))
meta = meta[meta$timepoint=="first_A",1:2]
meta$sampleID = gsub("_1", "", meta$sampleID)
table = merge(meta, table, by="sampleID")
urbanRural = factor(table$ruralUrban)
reads = log10(table$proportionReadsMapped)
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

par(mar=plot.mar)
boxplot(reads~urbanRural, 
        ylab="log10(proportion reads mapped)",# to CARD)",
        ylim = range(reads),
        outline=F,
        xaxt="n",
        cex.lab=2.5)
axis(side=1, at=c(1,2), labels=c("rural", "urban"), line=0, cex.axis=2.3)
x = ifelse(urbanRural=="rural", 1, 2)
points(reads~jitter(x), pch=16, col=ifelse(urbanRural=="rural", "blue", "red"), cex=2.2)
pval = t.test(reads~urbanRural)$p.value
print(pval) #0.01035906
mtext(getStars(pval), side=3, line=-5, cex=4.7)
mtext("D", side = 2, line = 2.58, cex=2.5, padj=0, las=1, at=-2.105) #logged
mtext("antibiotic\nresistance", side=3, line=0, cex=2.1)

###########
####virulence plot
table = read.table("MvirDB_propReads.txt", header=T, sep="\t", quote="",
                   colClasses=c("character", rep("numeric", 4)))
names(table)[1] = "sampleID"
##add metadata (same as antibiotic table antibiotic)
table = merge(meta, table, by="sampleID")
urbanRural = factor(table$ruralUrban)
reads = log10(table$proportionReadsMapped)

par(mar=plot.mar)
boxplot(reads~urbanRural, 
        ylab="log10(proportion reads mapped)",# to MvirDB)",
        ylim = range(reads),
        outline=F,
        xaxt="n",
        cex.lab=2.5)
axis(side=1, at=c(1,2), labels=c("rural", "urban"), line=0, cex.axis=2.3)
x = ifelse(urbanRural=="rural", 1, 2)
points(reads~jitter(x), pch=16, col=ifelse(urbanRural=="rural", "blue", "red"), cex=2.2)
pval = t.test(reads~urbanRural)$p.value
print(pval) #0.006257569
mtext(getStars(pval), side=3, line=-5, cex=4.7)
mtext("E", side = 2, line = 2.6, cex=2.5, padj=0, las=1, at=-1.26) #logged
mtext("virulence", side=3, line=1.2, cex=2.1)

dev.off()
