##make figure of metabolon data: PCoA, ROC and heatmap

rm(list=ls())
setwd("")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(vegan)
library(ROCR)
library(grid)

tiff("metabolonFig.tif", res=350, height=4666, width=5833)
layout(matrix(c(1,2,4,3,3,4), nrow=2, ncol=3, byrow=2), heights=c(1,1,2,1), 
       widths=c(1, .2, 1.2, 2))
plot.mar = c(4.2, 7, 4,.1) #plot margins
leg.mar = c(4.2, .1, 4,.1) #legend margins
par(cex.axis=1.5, cex.lab=1.5)

###########
###pcoa
##function to get get colors
getColors <- function(table) {
  colors = rep("blue", nrow(table)) #rural
  colors[table$ruralUrban=="urban"] = "red"
  return(colors)
}

table = read.table("MetabolonScaledImpData-transpose.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 337)))
names(table)[2] = "ruralUrban"
colors = getColors(table)
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
     col=colors, pch=16, cex=4, cex.lab=3)
mtext("A", side=3, line=1.4, cex=2.5, adj=0)

##write axes so can run model and get p values in legend
combine = cbind(table[,1:2], pcoa$CA$u)
write.table(combine, sep="\t", file="pcoaCorrected_metabolon.txt", quote=F, row.names=F, col.names=T)


###############
###pcoa legend
par(mar=leg.mar)
plot(1, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
legend("topleft", c("rural", "urban"), 
       col=c("blue", "red"), pch=16, cex=2.55)

##############
###ROC
levels = c("otu", "metabolon", "metadata", "metabolonOTU", "metadataOTU", "metabolonMetdata", "metabolonMetadataOTU")
colors = c("darkgreen", "purple", "orange", "black", "blue", "turquoise", "red")

##set up plot
par(mar=plot.mar)
plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate", ylab="True positive rate", cex.lab=2.5)
mtext("B", side=3, line=1.4, cex=2.5, adj=0)

##predictions
auc = rep(NA, length(levels))
for(i in 1:length(levels)) {
  anal = levels[i]
  table = read.table(paste("cforest_", anal, ".txt", sep=""), sep="\t", header=T, colClasses="numeric")
  pred = prediction(table$predicted, table$actual)
  perf = performance(pred, measure="tpr", x.measure="fpr")
  lines(x=unlist(slot(perf, "x.values")), y=unlist(slot(perf, "y.values")), col=colors[i], lwd=5)
  a = performance(pred,"auc")
  auc[i] = unlist(slot(a, "y.values"))
}

##legend
legend("bottomright", c("Input Data", "OTU", "Metabolite", "Metadata", "Metabolite+OTU", "Metadata+OTU", "Metabolite+Metadata", "Metabolite+OTU+Metadata          "),
       col=c("white", colors), lty=1, lwd=5, cex=2.4)
text(x=.97, y=seq(0.45, 0.02, length.out=length(auc)+1), labels=c("AUC", format(auc, digits=2)), cex=2.4)


#############
###heatmap
metabolon = read.table("metabolonModel_pValues.txt", sep="\t", header=T, colClasses=c("character", rep("numeric", 9)))
pathway = read.table("metabolonPathwayInfo.txt", sep="\t", header=T, colClasses="character")
pathway$biochemical = sub("*", "", pathway$biochemical, fixed=T) #* indicates compounds not officially confirmed
##adjust chemical names to have proper punctuation and pathway info
getPathway <- function(cheml) {
  if(cheml %in% pathway$biochemical) { #no need to modify
    r = which(pathway$biochemical==cheml)
    return(pathway$biochemical[r])
  } else {
    cheml = sub("^X", "", cheml) #need to remove X from chemicals that start with a number
    sp = strsplit(cheml, ".", fixed=T)[[1]] #split on dot, which replaces spaces, dashes, etc
    #find row that matches the split string
    r = grepl(sp[1], pathway$biochemical, fixed=T)
    i=2
    while(length(which(r)) > 1 & i <=length(sp)) {
      r2 = grepl(sp[i], pathway$biochemical, fixed=T)
      r = r & r2
      i = i+1
    }
    if(i > length(sp) & length(which(r)) > 1) { ##match not found
      cheml2 = sub(".", " ", cheml, fixed=T) #androsterone.sulfate has three hits - look for exact match (the other 2 have longer names that match)
      r = cheml2 == pathway$biochemical
      if(length(which(r))==1) {
        return(pathway$biochemical[r])
      } else {
        warning(paste("Multiple hits for chemical:", cheml))
        return(cheml)
      }
    } else {
      return(pathway$biochemical[r])
    }
  }
}

##sorted by pUrbanRural
metabolon = metabolon[order(metabolon$pValuesUrbanRural),]
sig = metabolon$adjustedPurbanRural < 0.05
###heatmap of average values
df = data.frame(names=metabolon$names[sig], meanRural=metabolon$meanRural[sig],meanUrban=metabolon$meanUrban[sig])
df$names = as.character(df$names)
for(i in 1:nrow(df)) {
  df$names[i] = getPathway(df$names[i])
}
df.m <- melt(df, id.vars="names") #columns to rows differentiated by measurement and value

gg <- ggplot(df.m, aes(x=variable, y=names)) + 
  geom_tile(aes(fill=value)) +
  scale_fill_gradient(low="white", high="purple4", limits=c(0, 2.2), name="mean\nscaled\nabundance") + 
  scale_x_discrete(name="", labels=c("rural", "urban")) +
  scale_y_discrete(name="", limits=rev(df$names)) +
  ggtitle("Metabolites") +
  theme_gray() +
  theme(axis.text.y=element_text(hjust=0, size=23, colour="black"), 
        axis.text.x=element_text(size=26, colour=c("blue", "red")),
        legend.title=element_text(size=18), legend.text=element_text(size=22), legend.key.height=unit(30, "points"), legend.key.width=unit(20, "points"),
        plot.title=element_text(size=20))

vp.top = viewport(layout=grid.layout(1, 2))
vp = viewport(layout.pos.col=2)
vps = vpTree(vp.top, vpList(vp))
plot(gg, vp=vps)

##new plot overlaid to get C in correct place
par(mar=c(.1,1,4,.1))
plot(1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
mtext("C", side=3, line=1.4, cex=2.5, adj=0)

dev.off()
