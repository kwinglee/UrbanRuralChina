##find correlations between OTU, metadata and metabolon data

rm(list=ls())
library(Kendall)

setwd("")

metadata = read.table("metadata_cleaned11-2015.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 55)))
metabolon = read.table("MetabolonScaledImpData-transpose.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 337)))

##get sample ids to be the same -> convert to numeric
names(metabolon)[1] = "sampleID"
names(metabolon)[2] = "ruralUrban"

categoricalVars = c("F9", "F13", "F17", "sex") #categorical metadata variables

##function that runs the comparisons
runComparisons <- function(table, splitCol, name) {
  namesOTU <- vector()
  namesOther <- vector()
  r <- vector()
  r.squared <- vector()
  pValues <- vector()
  pNonparam <- vector()
  nonparamTest <- vector()
  
  index <- 1
  
  pdf(paste("model_", name, "_plots.pdf",sep=""))
  
  for(i in 6:splitCol) {
    if(sum(table[,i] >0 ) > nrow(table) /4 && !all(table[,i]==1)) { 
     for ( j in (splitCol+1):ncol(table)) {
        if(!all(table[,j]==1)) { 
          namesOTU[index] <- names(table)[i]
          namesOther[index] <- names(table)[j]
          
          if(namesOTU[index] %in% categoricalVars | namesOther[index] %in% categoricalVars) {
            valI = table[,i]
            valJ = table[,j]
            if(namesOTU[index] %in% categoricalVars) {
              valI = table[,j]
              valJ = table[,i]
            }
            valJ = factor(valJ)
            
            r[index] = NA
            pValues[index] = t.test(valI~valJ)$p.value
            r.squared[index] = summary(lm(valI~valJ))$r.squared
            
            pNonparam[index] = wilcox.test(valI~valJ)$p.value
            nonparamTest[index] = "Wilcox"
            
            if(!is.na(pNonparam[index]) & pNonparam[index] < 0.05) { 
              myText <- paste(namesOTU[index], " vs. " ,namesOther[index] ,"\n", "p=" ,  
                              format(pValues[index] , digits=3), "r=", 
                              format(r[index], digits=3), "kendall p=" , 
                              format(pNonparam[index], digits=3)) 
              boxplot(valI~valJ, main=myText, xlab=namesOther[index], ylab=namesOTU[index])
              points(valI~valJ, pch=16, col = ifelse(table$ruralUrban == "rural", "blue", "red"))
            }
            
          } else {
            r[index] <- cor(table[,i], table[,j], method="spearman")
            aLm <- lm(table[,i] ~ table[,j])
            pValues[index] <- anova(aLm)$"Pr(>F)"[1]
            r.squared[index] = summary(aLm)$r.squared
            pNonparam[index] <- Kendall(table[,i], table[,j])$sl[1]
            nonparamTest[index] <- "Kendall"
            
            if((!is.na(pNonparam[index]) & pNonparam[index] < 0.05) |
               (pValues[index] < 0.05)) { 
              myText <- paste(namesOTU[index], " vs. " ,namesOther[index] ,"\n", "p=" ,  
                              format(pValues[index] , digits=3), "r=", 
                              format(r[index], digits=3), "kendall p=" , 
                              format(pNonparam[index], digits=3)) 
              plot(table[,j],table[,i] , main=myText, pch=16, xlab=namesOther[index], ylab=namesOTU[index],
                   col = ifelse(table$ruralUrban == "rural", "blue", "red"))
              abline(aLm)
            }
          }
          
          index <- index + 1
        }
      }
    }
  }
  
  dev.off()
  dFrame <- data.frame( namesOTU, namesOther, pLM=pValues, r, r.squared, nonparamTest, pNonparam) 
  dFrame <- dFrame [order(dFrame$pLM),]
  dFrame$adjLMp <-  p.adjust(dFrame$pLM, method = "BH")
  dFrame$adjNonparmP <- p.adjust(dFrame$pNonparam, method = "BH")
  
  write.table( file= paste("model_", name, "_pValues.txt", sep=""), dFrame, row.names=FALSE, sep="\t")
}


taxaLevels <- c("phylum","class","order","family","genus", "otu")


for(taxa in taxaLevels ) {
  print(taxa)
  inFileName <- paste(taxa,  "_taxaAsColumnsLogNorm_WithMetadata.txt", sep ="")
  if(taxa=="otu") {
    inFileName <- "abundantOTUForwardTaxaAsColumnsLogNormalWithMetadata.txt"
  }
  myT <-read.table(inFileName,header=TRUE,sep="\t")
  numCols <- ncol(myT)
  myColClasses <- c(rep("character",5), rep("numeric", numCols-5))
  if(taxa=="otu") { #missing read number
    myColClasses <- c(rep("character",4), rep("numeric", numCols-4))
  }
  myT <-read.table(inFileName,header=TRUE,sep="\t",colClasses=myColClasses)
  
  ##add read number to OTU
  if(taxa=="otu") {
    readNumber = rep("1", nrow(myT))
    myT <- cbind(myT[,1], readNumber, myT[,-1])
    names(myT)[1] = "sampleID"
    numCols = ncol(myT)
  }
  
  ## only the forward reads
  myT <- myT[ which( myT$readNumber == "1"), ]
  
  ##metadata and metabolite are only first timepoint
  myT = myT[myT$timepoint == "first_A",]
  
  ##convert sample IDs to numeric like metadata and metabolon
  names = myT$sampleID
  names = sub("^0", "", names)
  names = sub("A_1$", "", names)
  names = sub("B", "", names)
  names = sub("A", "", names)
  myT$sampleID = as.numeric(names)
  
  ##metadata vs OTU
  mrg.md = merge(myT, metadata, by=c("sampleID", "ruralUrban"))
  runComparisons(mrg.md, numCols, paste("metadata_v_", taxa, sep=""))
  
  ##metabolite dataset
  mrg.mb = merge(myT, metabolon, by=c("sampleID", "ruralUrban"))
  runComparisons(mrg.mb, numCols, paste("metabolon_v_", taxa, sep=""))
}

##metadata vs. metabolon
##add three extra columns to beginning so can feed it through same method
mrg = merge(metabolon, metadata, by=c("sampleID", "ruralUrban"))
X = rep(NA, nrow(mrg))
mrg = cbind(X, X, X, mrg)
runComparisons(mrg, ncol(metabolon)+3, "metadata_v_metabolon")
