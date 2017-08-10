##make supplemental figure comparing timepoints (16S)
##for genera and OTU, draw time 1 vs. time 2 for MDS1, MDS2

rm(list=ls())

setwd("")

genmar = c(4, 5.5, 3, .1) #genus margin top large for title
otumar = c(5.5, 5.5, 1.5, .1) #otu bottom margin large for legend
tiff("T1vT2SuppFig.tif", res=300, height=2800, width=2800)
par(mfrow=c(2,2), cex=1)

##draw graph colored by urban-rural, x is time 1 values, y is time 2
drawGraph <- function(rurUrb, x, y) {
  lm=lm(y~x)
  r = summary(lm)$r.squared
  plot(x=x, y=y, xlab="first timepoint", ylab="second timepoint", main="", 
       col=ifelse(rurUrb == "rural", "blue", "red"), pch=15)
  abline(lm)
  legend("bottomright", format(r, digits=2), bty="n", cex=1.5)
}

labelcex = 1.5

##plot MDS1 and MDS2 time 1 vs time 2 for the given taxa
mdsGraph <- function(taxa) {
  fileName <- paste("pcoaCorrected_", taxa, ".txt", sep ="")
  table <-read.table(fileName,header=TRUE,sep="\t")
  numCols <- ncol(table)
  cc=c("character", "numeric", "numeric", "character", "character", rep("numeric", numCols-5))
  table <-read.table(fileName,header=TRUE,sep="\t",colClasses=cc)
  table = table[table$readNumber==1,]
  
  tab1 = table[table$timepoint=="first_A",]
  tab2 = table[table$timepoint=="second_B",]
  df1 = data.frame(patientID=tab1$patientID, ruralUrban=tab1$ruralUrban, MDS1A=tab1$MDS1, MDS2A=tab1$MDS2)
  df2 = data.frame(patientID=tab2$patientID, ruralUrban=tab2$ruralUrban, MDS1B=tab2$MDS1, MDS2B=tab2$MDS2)
  mrg = merge(df1, df2, by=c("patientID", "ruralUrban"))
  
  if(taxa=="genus") {
    par(mar=genmar)
  } else {
    par(mar=otumar)
  }
  
  drawGraph(mrg$ruralUrban, mrg$MDS1A, mrg$MDS1B)
  mtext(taxa, side=2, line=4, cex=labelcex)
  if(taxa == "genus") {
    mtext("MDS1", side=3, line=1, cex=labelcex)
  }
  drawGraph(mrg$ruralUrban, mrg$MDS1A, mrg$MDS1B)
  if(taxa == "genus") {
    mtext("MDS2", side=3, line=1, cex=labelcex)
  }
}

taxaLevels <- c("genus", "OTU")
for(i in 1:length(taxaLevels)) {
  taxa = taxaLevels[i]
  mdsGraph(taxa)
}

##legend
par(oma=c(0,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(0,0,type="n", bty="n", xaxt="n", yaxt="n")
legend("bottom", inset=c(0,0), horiz=T,
       legend=c("rural", "urban"), 
       col = c("blue" ,"red"),
       pch = 15, cex=1.2)
dev.off()