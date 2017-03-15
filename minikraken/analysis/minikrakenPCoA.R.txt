##PCoA of minikraken analysis

rm(list=ls())

library(vegan)

setwd("")

levels = c("domain", "phylum", "class", "order", "family", "genus", "species")

pdf("minikrakenPCoA.pdf")
for(lev in levels) {
  print(lev)
  file = paste("minikraken_merged_taxaAsCol_logNorm_", lev, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  pcoa = capscale(table[,-(1:3)]~1,distance="bray")
  
  ##get percent variance
  eig = eigenvals(pcoa)
  if(any(eig<0)) {#need to check there are no negatives
    warning(paste("NEGATIVE EIGENVALUES-percent variance is incorrect for ", lev, sep=""))
  }
  var = eig/sum(eig)*100 
  
  colors = rep("blue", nrow(table)) #rural
  colors[table$ruralUrban=="urban"] = "red"
  
  par(mar=c(4, 4, 4, 6), xpd=T)
  plot(x=pcoa$CA$u[,1], y=pcoa$CA$u[,2], 
       xlab=paste("MDS1 (", format(var[[1]], digits=0), "%)", sep=""), 
       ylab=paste("MDS2 (", format(var[[2]], digits=0), "%)", sep=""), 
       main=lev, col=colors, pch=16)
  legend("topright", inset=c(-.2, 0),
         legend=c("rural", "urban"),
         col=c("blue", "red"), pch=16)
  
  ##write axes
  combine = cbind(table[,1:3], pcoa$CA$u)
  write.table(combine, sep="\t", file=paste("minikraken_PcoaCorrected_", lev, ".txt",sep=""), quote=F, row.names=F, col.names=T)
}
dev.off()