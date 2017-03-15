##HUMAnN volcano plots with p-values

rm(list=ls())

setwd("")

levels = c("module", "pathway")

tiff("humann_volcanoPlot.tif", height=1300, width=2500, res=300)
par(oma=c(3,0,0,0), mfrow=c(1,2), mar = c(5.1, 4.5, 4.1, 2.1))

for(i in 1:length(levels)) {
  lev = levels[i]
  print(lev)
  table = read.table(paste("humann_otuModel_pValues_unlog_", lev, ".txt",sep=""), sep="\t", header=T)
  
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
       main = paste(lev, "\np =", format(p, digits = 2)), cex.lab=1.5, cex.main=1.5, col=colors, cex=.5)
  abline(h=-log10(0.05), lty=3, col="gray")
  abline(v=0, lty=3, col="gray")
}

##legend
par(oma=c(0,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", cex=1, ncol=2,
       legend=c("taxa significantly more abundant in rural", "taxa more abundant in rural",
                "taxa significantly more abundant in urban", "taxa more abundant in urban"),
       col=c("blue", "lightblue", "red", "pink"),
       pch=16)
dev.off()