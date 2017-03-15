##metadata PCoA

rm(list=ls())

library("vegan")

setwd("")

tiff("metadataFig.tiff", res=300, height=1500, width=1800)
par(cex.axis=1.5, cex.lab=1.5)

###pcoa
##function to get get colors
getColors <- function(table) {
  colors = rep("blue", nrow(table)) #rural
  colors[table$ruralUrban=="urban"] = "red"
  return(colors)
}

table = read.table("metadata_cleaned11-2015.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 55)))
colors = getColors(table)
pcoa <- capscale(table[,-(1:2)]~1,distance="euclidean")

##get percent variance
eig = eigenvals(pcoa)
if(any(eig<0)) {#need to check there are no negatives
  warning(paste("NEGATIVE EIGENVALUES-percent variance is incorrect for ", taxa, sep=""))
}
var = eig/sum(eig)*100

par(xpd=T, mar=c(4, 4.2, 2, 6.5))
plot(x=pcoa$CA$u[,1], y=pcoa$CA$u[,2],  
     xlab=paste("MDS1 (", format(var[[1]], digits=0), "%)", sep=""), 
     ylab=paste("MDS2 (", format(var[[2]], digits=0), "%)", sep=""), 
     col=colors, pch=16, cex=2)
  
legend("topright", c("rural", "urban"), inset=c(-.31,0),
       col=c("blue", "red"), pch=16, cex=1.5)

##write axes so can run model and get p values in legend
combine = cbind(table[,1:2], pcoa$CA$u)
write.table(combine, sep="\t", file="pcoaCorrected_metadata.txt", quote=F, row.names=F, col.names=T)

dev.off()
