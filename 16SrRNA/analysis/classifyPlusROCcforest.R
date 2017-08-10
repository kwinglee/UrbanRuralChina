##run classifiers on 16S data using cforest

rm(list=ls())
setwd("")

library(party)
library(ROCR)

##function to run permutations from given predictions and true class, adding the lines to the plot
permut <- function(predClass, trueClass) {
  for(i in 1:30) {
    pred = prediction(sample(predClass, size=length(predClass)), trueClass)
    perf = performance(pred, measure="tpr", x.measure="fpr")
    lines(x=unlist(slot(perf, "x.values")), y=unlist(slot(perf, "y.values")), col="gray")
  }
}

##function that classifies the input table using bagging then draws ROC curve in the color col and runs permutations and returns the AUC
##and writes the predictions to file with name so don't have to keep running classifier
classifyROC <- function(table, col, name) {
  ##convert rural/urban to -1/1
  names(table)[1] = "ruralUrban"
  class = rep(-1, nrow(table)) #rural
  class[table$ruralUrban == "urban"] = 1 #urban
  
  dat = cbind(class, table[,-1])
  
  ##get predictions
  bag.pred = rep(NA, nrow(table))
  for(i in 1:nrow(table)) {
    bag = cforest(class ~ ., data=dat[-i,])
    bag.pred[i] = predict(bag, newdata=dat[i,-1])
  }
  
  ##write results
  write.table(data.frame(actual=class, predicted=bag.pred), paste("cforest_", name, ".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
  
  ##run permutations
  # permut(bag.pred, class)
  
  ##draw ROC curve
  pred.rocr = prediction(bag.pred, class)
  perf.rocr = performance(pred.rocr, measure="tpr", x.measure="fpr")
  lines(x=unlist(slot(perf.rocr, "x.values")), y=unlist(slot(perf.rocr, "y.values")), col=col, lwd=3)
  
  ##return AUC
  auc = performance(pred.rocr,"auc")
  return(unlist(slot(auc, "y.values")))
}

tiff("ROCcforest-microbiome.tiff", res=300, height=1500, width=2000)

par(xpd=TRUE, mar=c(4.1, 4.1, 1.1, 0.1)) #margins are bottom, left, top, right; default is 5.1, 4.1, 4.1, 2.1
plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate", ylab="True positive rate")

##otu
otu.auc = classifyROC(read.table("abundantOTULogNormal.txt", header=T, sep="\t", row.names=1, colClasses=c("character", "character", rep("numeric", 867))), 
                      "blue", "otu")
##genera
gen.auc = classifyROC(read.table("genus_forSVMLight.txt", header=T, sep="\t", row.names=1, colClasses=c("character", "character", rep("numeric", 347))), 
                      "green", "genus")
##family
fam.auc = classifyROC(read.table("family_forSVMLight.txt", header=T, sep="\t", row.names=1, colClasses=c("character", "character", rep("numeric", 126))), 
                      "red", "family")
##order
ord.auc = classifyROC(read.table("order_forSVMLight.txt", header=T, sep="\t", row.names=1, colClasses=c("character", "character", rep("numeric", 56))), 
                      "purple", "order")
##class
cla.auc = classifyROC(read.table("class_forSVMLight.txt", header=T, sep="\t", row.names=1, colClasses=c("character", "character", rep("numeric", 38))), 
                      "turquoise", "class")
##phylum
phy.auc = classifyROC(read.table("phylum_forSVMLight.txt", header=T, sep="\t", row.names=1, colClasses=c("character", "character", rep("numeric", 21))), 
                      "black", "phylum")

legend("bottomright",
       legend = c(paste("Input Data", "AUC", sep="  "), paste("OTU", format(otu.auc, digits=2), sep="          "), paste("Genera", format(gen.auc, digits=2), sep="      "), 
                  paste("Family", format(fam.auc, digits=2), sep="       "), paste("Order", format(ord.auc, digits=2), sep="         "), paste("Class", format(cla.auc, digits=2), sep="         "),
                  paste("Phyla", format(phy.auc, digits=2), sep="         ")),#, "permutations"), 
       col=c("white", "blue", "green", "red", "purple", "turquoise", "black"), lty=1, lwd=2, cex=1)

dev.off()

tiff("ROCcforest-otuMetaboliteMetadata.tiff", res=300, height=1500, width=2000)

par(xpd=TRUE, mar=c(4.1, 4.1, 1.1, 0.1)) #margins are bottom, left, top, right; default is 5.1, 4.1, 4.1, 2.1
plot(1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate", ylab="True positive rate")

##otu
otu.auc = classifyROC(read.table("abundantOTULogNormal.txt", header=T, sep="\t", row.names=1, colClasses=c("character", "character", rep("numeric", 867))), 
                      "blue", "otu")
##metabolon
mb.auc = classifyROC(read.table("MetabolonScaledImpData-transpose.txt", header=T, sep="\t", row.names=1, colClasses=c("numeric", "character", rep("numeric", 337))), 
                     "black", "metabolon")
##metadata
md.auc = classifyROC(read.table("metadata_cleaned11-2015.txt", header=T, sep="\t", row.names=1, colClasses=c("numeric", "character", rep("numeric", 55))), 
                     "red", "metadata")
##metabolon+OTU
mbotu.auc = classifyROC(read.table("MetabolonWithOTU.txt", header=T, sep="\t", row.names=1, colClasses=c("numeric", "character", rep("numeric", 1204))), 
                        "green", "metabolonOTU")
##metadata+OTU
mdotu.auc = classifyROC(read.table("OTUWithMetadata_v3.txt", header=T, sep="\t", row.names=1, colClasses=c("numeric", "character", rep("numeric", 922))), 
                        "purple", "metadataOTU")
##metabolon+metadata
mdmb.auc = classifyROC(read.table("MetabolonWithMetadata_v3.txt", header=T, sep="\t", row.names=1, colClasses=c("numeric", "character", rep("numeric", 392))), 
                       "turquoise", "metabolonMetdata")
##metabolon+metadata+OTU
all.auc = classifyROC(read.table("OTUMetadataMetabolon_v3.txt", header=T, sep="\t", row.names=1, colClasses=c("numeric", "character", rep("numeric", 1259))), 
                      "orange", "metabolonMetadataOTU")


legend("bottomright", c(paste("Input Data", "AUC", sep="                          "), paste("Metabolite", format(mb.auc, digits=2), sep="                          "), paste("OTU", format(otu.auc, digits=2), sep="                           "), 
                        paste("Metadata", format(md.auc, digits=2), sep="                             "), paste("Metabolite+OTU", format(mbotu.auc, digits=2), sep="                  "), paste("Metadata+OTU", format(mdotu.auc, digits=2), sep="                    "),
                        paste("Metabolite+Metadata", format(mdmb.auc, digits=2), sep="          "), paste("Metabolite+OTU+Metadata", format(all.auc, digits=2), sep=" ")),# "permutations"), 
       col=c("white", "black", "blue", "red", "green", "purple", "turquoise", "orange"), lty=1, lwd=2, cex=0.75)

dev.off()