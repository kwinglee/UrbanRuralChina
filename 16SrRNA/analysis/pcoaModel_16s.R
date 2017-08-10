##mixed linear models of MDS ~ time + urbanRural + 1|subjectID for 16S PCoA axes

rm(list=ls())

library("pscl")
library("lmtest")
library("nlme")

setwd("")

taxaLevels <- c("phylum","class","order","family","genus", "otu")

for(taxa in taxaLevels ) {
  print(taxa)
  inFileName <- paste("pcoaCorrected_", taxa,  ".txt", sep ="")
  myT <-read.table(inFileName,header=TRUE,sep="\t")
  numCols <- ncol(myT)
  myColClasses <- c("character", "numeric", "numeric", "character", "character", rep("numeric", numCols-5))
  myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)
  
  names <- vector()
  pValuesTime <- vector()
  pValuesSubject <- vector()
  pValuesUrbanRural <- vector()
  meanMDS <- vector()
  meanUrban1 <- vector()
  meanUrban2 <- vector()
  meanRural1 <- vector()
  meanRural2 <- vector()
  pValuesUrbanRuralWilcoxT1 <- vector()
  pValuesUrbanRuralWilcoxT2 <- vector()
  index <- 1
  pdf( paste("pcoaModel_boxplots_", taxa, ".pdf", sep=""))
  
  for( i in 6:numCols) {
      
      mds <- myT[,i]
      meanMDS[index] <- mean(mds)
      meanUrban1[index] <- mean(mds[myT$ruralUrban=="urban" & myT$timepoint=="first_A"])
      meanRural1[index] <- mean(mds[myT$ruralUrban=="rural" & myT$timepoint=="first_A"])
      meanUrban2[index] <- mean(mds[myT$ruralUrban=="urban" & myT$timepoint=="second_B"])
      meanRural2[index] <- mean(mds[myT$ruralUrban=="rural" & myT$timepoint=="second_B"])
      time <- factor(myT$timepoint)
      patientID <- factor(myT$patientID )
      urbanRural <- factor(myT$ruralUrban)
      
      myFrame <- data.frame(mds, time, patientID, urbanRural )
      
      fullModel <- gls( mds~  time + urbanRural , method="REML",correlation=corCompSymm(form=~1|factor(patientID)),
                        data = myFrame )
      reducedModel <- gls( mds~  time + urbanRural , method="REML",	data = myFrame )
      fullModelLME <- lme(mds~  time + urbanRural , method="REML", random = ~1|factor(patientID), data = myFrame)		
      
      pValuesTime[index] <- anova(fullModelLME)$"p-value"[2]
      pValuesUrbanRural[index] <- anova(fullModelLME)$"p-value"[3]
      pValuesSubject[index] <-  anova(fullModelLME, reducedModel)$"p-value"[2]
      intraclassCoefficient<- coef(fullModel$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]
      names[index] = names(myT)[i]
      
      graphMain =  paste( names(myT)[i], " pTime=", format(pValuesTime[index], digits=3), "\n",
                          " pRuralUrban= ", format(pValuesUrbanRural[index],digits=3), 
                          " pSubject= " , format(	pValuesSubject[index], digits=3), "\n",
                          " icc= " , format( intraclassCoefficient, digits=3 ), sep="")
      
      par(mfrow=c(3,1))
      
      plot( mds[myT$timepoint=="first_A"] ~ urbanRural[myT$timepoint=="first_A"], ylab = names[index],
            main = graphMain )	
      
      stripchart(mds[myT$timepoint=="first_A"] ~ urbanRural[myT$timepoint=="first_A"], 
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE )		
      
      plot( mds[myT$timepoint=="second_B"] ~ urbanRural[myT$timepoint=="second_B"], ylab=names[index])	
      
      stripchart(mds[myT$timepoint=="second_B"] ~ urbanRural[myT$timepoint=="second_B"], 
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])		
      
      plot(mds ~ patientID , las=3 , ylab = names[index]) 		
      
      stripchart(mds[myT$ruralUrban=="rural"] ~ patientID[myT$ruralUrban=="rural"], 
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE,	col = "blue")
      
      stripchart(mds[myT$ruralUrban=="urban"] ~ patientID[myT$ruralUrban=="urban"], 
                 data = myFrame,vertical = TRUE, pch = 21, add=TRUE,	col = "red")
      
      ##non parametric test
      time1 = myT$timepoint=="first_A"
      bug1 = myT[time1, i]
      urb1 = factor(myT$ruralUrban[time1])
      time2 = myT$timepoint=="second_B"
      bug2 = myT[time2, i]
      urb2 = factor(myT$ruralUrban[time2])
      
      pValuesUrbanRuralWilcoxT1[index] = wilcox.test(bug1~urb1, exact=T)$p.value
      pValuesUrbanRuralWilcoxT2[index] = wilcox.test(bug2~urb2, exact=T)$p.value
      #note, if don't give exact=F, get warning: cannot compute exact p-value with ties
      
      index=index+1
      
    }
  
  dFrame <- data.frame( names,meanMDS,meanUrban1, meanUrban2, meanRural1, meanRural2, pValuesTime ,pValuesSubject, pValuesUrbanRural, pValuesUrbanRuralWilcoxT1, pValuesUrbanRuralWilcoxT2)
  dFrame$UrbanToRural <- ((meanUrban1 + meanUrban2)/2) / ((meanRural1 + meanRural2)/2)
  dFrame$adjustedPtime <- p.adjust( dFrame$pValuesTime, method = "BH" )
  dFrame$adjustedPsubject <- p.adjust( dFrame$pValuesSubject, method = "BH" )
  dFrame$adjustedPurbanRural <- p.adjust( dFrame$pValuesUrbanRural, method = "BH" )
  dFrame$adjustedPurbanRuralWilcoxT1 <- p.adjust(dFrame$pValuesUrbanRuralWilcoxT1, method="BH")
  dFrame$adjustedPurbanRuralWilcoxT2 <- p.adjust(dFrame$pValuesUrbanRuralWilcoxT2, method="BH")
  dFrame <- dFrame [order(dFrame$adjustedPurbanRural),]
  write.table(dFrame, file=paste("pcoaModel_pValues_", taxa, ".txt",sep=""), sep="\t",row.names=FALSE)
  dev.off()
}