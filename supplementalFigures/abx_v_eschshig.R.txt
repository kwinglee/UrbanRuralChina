##16S and minikraken Proteobacteria-Escherichia/Shigella vs. proportion antibiotic resistance reads

rm(list=ls())
setwd("")

##get antibiotics numbers
table = read.table("pro_homolog_results.txt", header=T, sep="\t", stringsAsFactors = F, quote="")
abx = data.frame(sampleID=table$sample, prop=table$proportionReadsMapped, stringsAsFactors = F)

tiff("abx_v_eschshig.tif", height=2000, width=6000, res=300)
par(mfrow=c(2,6))

plot.mar = c(5.1, 4.2, 4.1, 1.5)

####16S
##function that takes the given taxonomic level and name of the taxa of interest and plots it with the given
##label lab in the top left
plot16S <- function(level, taxa, lab="") {
  fileName = paste("level,  "_taxaAsColumnsLogNorm_WithMetadata.txt", sep ="")
  table = read.table(fileName,header=TRUE,sep="\t")
  numCols = ncol(table)
  ColClasses = c(rep("character",5), rep("numeric", numCols-5))
  table = read.table(fileName,header=TRUE,sep="\t",colClasses=ColClasses)
  
  ##forward read only
  ##first time point only
  ##fix sample names
  table = table[table$readNumber=="1" & table$timepoint=="first_A",]
  table$sampleID = gsub("_1", "", table$sampleID)
  
  ##merge
  mrg = merge(abx, table, by="sampleID")
  abun = mrg[,names(mrg)==taxa]
  
  ##get r and p
  pVal = read.table(paste("cardsAssoc_ProHomolog_propMapped_v_16S_", level, "_pValues.txt", sep=""),
                    header=T, sep="\t", stringsAsFactors = F)
  p = pVal$adjLMp[pVal$namesOther==taxa]
  r = pVal$r[pVal$namesOther==taxa]
  
  taxa = gsub(".", "/", taxa, fixed=T) #fixed genus
  
  ##plot
  par(mar=plot.mar)
  col = ifelse(mrg$ruralUrban=="rural", "blue", "red")
  plot(x=abun, y=mrg$prop,
       xlab = "16S log relative abundance",
       ylab = "proportion reads mapped to CARD",
       col = col,
       pch = 16,
       cex=2, cex.lab=1.3)
  if(level=="genus") {
    title(taxa, cex.main=1.5, font.main=4)
  } else {
    title(taxa, cex.main=1.5)
  }
  mtext(paste("(", level, ")", sep=""), side=3, cex=1.1)
  legend("topleft",
         legend=c(paste("r=", format(r, digits=2), sep=""),
                  paste("p=", format(p, digits=2), sep="")), cex=1.2)
  abline(lm(mrg$prop~abun))
  mtext(lab, side=3, line=1.5, cex=2, adj=0)
}

plot16S("phylum", "Proteobacteria", "A")
plot16S("class", "Gammaproteobacteria")
plot16S("order", "Enterobacteriales")
plot16S("family", "Enterobacteriaceae")
plot16S("genus", "Escherichia.Shigella")

###blank plot with key
par(mar=c(.1, 4.1, 4.1, .1))
plot(1, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
legend("topleft",
       title = expression(underline("Key")),
       legend=c("rural", "urban", "linear regression line"),
       col=c("blue", "red", "black"),
       pch=c(16, 16, NA), 
       lty=c(NA, NA, 1),
       cex=1.5)

####minikraken
##function that takes the given taxonomic level and name of the taxa of interest and plots it with the given
##label lab in the top left
plotMinikraken <- function(level, taxa, lab="") {
  file = paste("minikraken_merged_taxaAsCol_logNorm_", level, ".txt", sep="")
  table = read.table(file, header=T, sep="\t")
  ncol = ncol(table)
  table = read.table(file, header=T, sep="\t", 
                     colClasses=c("character", "numeric", "character", rep("numeric", ncol-3)))
  
  ##merge
  mrg = merge(abx, table, by="sampleID")
  abun = mrg[,names(mrg)==taxa]
  
  ##get r and p
  pVal = read.table(paste("cardsAssoc_ProHomolog_propMapped_v_minikraken_", level, "_pValues.txt", sep=""),
                    header=T, sep="\t", stringsAsFactors = F)
  p = pVal$adjLMp[pVal$namesOther==taxa]
  r = pVal$r[pVal$namesOther==taxa]
  
  ##plot
  par(mar=plot.mar)
  col = ifelse(mrg$ruralUrban=="rural", "blue", "red")
  plot(x=abun, y=mrg$prop,
       xlab = "WGS log relative abundance",
       ylab = "proportion reads mapped to CARD",
       col = col,
       pch = 16,
       cex=2, cex.lab=1.3)
  if(level=="genus") {
    title(taxa, cex.main=1.5, font.main=4)
  } else {
    title(taxa, cex.main=1.5)
  }
  mtext(paste("(", level, ")", sep=""), side=3, cex=1.1)
  legend("topleft",
         legend=c(paste("r=", format(r, digits=2), sep=""),
                  paste("p=", format(p, digits=2), sep="")), cex=1.2)
  abline(lm(mrg$prop~abun))
  mtext(lab, side=3, line=1.5, cex=2, adj=0)
}

plotMinikraken("phylum", "Proteobacteria", "B")
plotMinikraken("class", "Gammaproteobacteria")
plotMinikraken("order", "Enterobacteriales")
plotMinikraken("family", "Enterobacteriaceae")
plotMinikraken("genus", "Escherichia")
plotMinikraken("genus", "Shigella")

dev.off()