##generate supplemental figure of ROC curves from classifications using the rpart function

rm(list=ls())
setwd("")

library(ROCR)

tiff("rpartROCSuppFig.tif", res=300, height=4000, width=4000)

par(cex.axis=1.5, cex.lab=1.5, xpd=T, mfrow=c(2,1), mar=c(4,4,3,25))

############
####microbiome
taxaLevels = c("phylum", "class", "order", "family", "genus", "otu")
colors = c("blue", "black", "red", "purple", "turquoise", "darkgreen")

##set up plot
plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate", ylab="True positive rate")
mtext("A", side=3, line=1.5, cex=2, adj=0)

##predictions
auc = rep(NA, length(taxaLevels))
for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  table = read.table(paste("rpart_", taxa, ".txt", sep=""), sep="\t", header=T, colClasses="numeric")
  pred = prediction(table$predicted, table$actual)
  perf = performance(pred, measure="tpr", x.measure="fpr")
  lines(x=unlist(slot(perf, "x.values")), y=unlist(slot(perf, "y.values")), col=colors[i], lwd=4)
  a = performance(pred,"auc")
  auc[i] = unlist(slot(a, "y.values"))
}

legend("topright", inset=c(-.375,0),
       legend = c("Input Data", "OTU", "Genus", "Family", "Order", "Class", "Phylum                "),
       col=c("white", rev(colors)), lty=1, lwd=2, cex=1.5, bty="n")
legend("topright", inset=c(-.63,0),
       legend=c("AUC", rev(format(auc, digits=2))), cex=1.5, bty="n")

###########
###metabolite/metadata
levels = c("otu", "metabolon", "metadata", "metabolonOTU", "metadataOTU", "metabolonMetadata", "metabolonMetadataOTU")
colors = c("darkgreen", "purple", "orange", "black", "blue", "turquoise", "red")

##set up plot
plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate", ylab="True positive rate")
mtext("B", side=3, line=1.5, cex=2, adj=0)

##predictions
auc = rep(NA, length(levels))
for(i in 1:length(levels)) {
  anal = levels[i]
  table = read.table(paste("rpart_", anal, ".txt", sep=""), sep="\t", header=T, colClasses="numeric")
  pred = prediction(table$predicted, table$actual)
  perf = performance(pred, measure="tpr", x.measure="fpr")
  lines(x=unlist(slot(perf, "x.values")), y=unlist(slot(perf, "y.values")), col=colors[i], lwd=4)
  a = performance(pred,"auc")
  auc[i] = unlist(slot(a, "y.values"))
}

##legend
legend("topright", inset=c(-.61,0),
       c("Input Data", "OTU", "Metabolite", "Metadata", "Metabolite+OTU", "Metadata+OTU", "Metabolite+Metadata", "Metabolite+OTU+Metadata          "),
       col=c("white", colors), lty=1, lwd=2, cex=1.5, bty="n")
legend("topright", inset=c(-.63,0), 
       legend=c("AUC", format(auc, digits=2)), cex=1.5, bty="n")

dev.off()
