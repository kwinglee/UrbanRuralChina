##PCoA of HUMAnN analysis (unlogged data)

rm(list=ls())
setwd("")

library(vegan)

levels = c("module", "pathway")

pdf("humann_pcoa_unlog.pdf")
for(lev in levels) {
  print(lev)
  file = paste("humann_keggAsCol_withRurUrb_", lev, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", colClasses=c(rep("character", 2), rep("numeric", ncol-2)))
  
  pcoa = capscale(table[,-(1:6)]~1,distance="bray")
  
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
  combine = cbind(table[,1:2], pcoa$CA$u)
  write.table(combine, sep="\t", file=paste("humann_PcoaCorrected_unlog_", lev, ".txt",sep=""), quote=F, row.names=F, col.names=T)
}
dev.off()