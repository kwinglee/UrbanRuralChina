##minikraken volcano plots with p-values

rm(list=ls())

setwd("")

taxaLevels = c("phylum", "class", "order", "family", "genus", "species")

tiff("minikraken_volcanoPlot.tif", height=2400, width=3000, res=300)
par(oma=c(4,0,0,0), mfrow=c(2,3))

for(i in 1:length(taxaLevels)) {
  lev = taxaLevels[i]
  print(lev)
  table = read.table(paste("minikraken_otuModel_pValues_", lev, ".txt",sep=""), sep="\t", header=T)
  
  ##contingency table
  ctbl = cbind(hiUrban = c(sum(table$adjustedPurbanRural <= 0.05 & table$UrbanToRural > 1),
                           sum(table$adjustedPurbanRural > 0.05 & table$UrbanToRural > 1)),
               hiRural = c(sum(table$adjustedPurbanRural <= 0.05 & table$UrbanToRural < 1),
                           sum(table$adjustedPurbanRural > 0.05 & table$UrbanToRural < 1)))
  row.names(ctbl) = c("sig", "notSig")
  print(ctbl)
  ##chi square test
  p = chisq.test(ctbl, simulate.p.value = T)$p.value
  
  ##volcano plot
  colors = rep("black", nrow(table))
  colors[table$meanUrban > table$meanRural & table$adjustedPurbanRural < 0.05] = "red"
  colors[table$meanUrban > table$meanRural & table$adjustedPurbanRural > 0.05] = "pink"
  colors[table$meanUrban < table$meanRural & table$adjustedPurbanRural < 0.05] = "blue"
  colors[table$meanUrban < table$meanRural & table$adjustedPurbanRural > 0.05] = "lightblue"
  plot(x = log2(table$UrbanToRural), y = -1*log10(table$adjustedPurbanRural),
       xlab = "log2 urban/rural fold change", ylab= "-log10 adjusted p-values", pch=16,
       main = paste(lev, "\np =", format(p, digits = 2)), cex.lab=1.5, cex.main=1.5, col=colors)
  abline(h=-log10(0.05), lty=3, col="gray")
  abline(v=0, lty=3, col="gray")
}

##legend
par(oma=c(0,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", cex=1.5, ncol=2,
       legend=c("taxa significantly more abundant in rural", "taxa more abundant in rural",
                "taxa significantly more abundant in urban", "taxa more abundant in urban"),
       col=c("blue", "lightblue", "red", "pink"),
       pch=16)
dev.off()