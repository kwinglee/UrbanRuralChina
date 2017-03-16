##model results for metabolite data using linear models and nonparametric (Wilcox) test

rm(list=ls())

setwd("")

metabolon = read.table("MetabolonScaledImpData-transpose.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 337)))

###p-values
##output vectors
names <- vector()
pValuesUrbanRural <- vector()
meanMetabolite <- vector()
sdMetabolite <- vector()
meanUrban <- vector()
sdUrban <- vector()
meanRural <- vector()
sdRural <- vector()
pWilcox <- vector()
r.squared <- vector()
index <- 1
urbanRural <- factor(metabolon$classification)
  
for( i in 3:ncol(metabolon)) {
  if(sum(metabolon[,i] != 0 ) > nrow(metabolon) / 4  && !all(metabolon[,i]==1)) { 
    met <- metabolon[,i]
    meanMetabolite[index] <- mean(met)
    meanUrban[index] <- mean(met[metabolon$classification=="urban"])
    meanRural[index] <- mean(met[metabolon$classification=="rural"])
    sdMetabolite[index] <- sd(met)
    sdUrban[index] <- sd(met[metabolon$classification=="urban"])
    sdRural[index] <- sd(met[metabolon$classification=="rural"])
    
    myFrame <- data.frame(met, urbanRural)
      
    model = lm(met ~ urbanRural)
    names[index] = names(metabolon)[i]
    pValuesUrbanRural[index] <- anova(model)$`Pr(>F)`[1]
    r.squared[index] = summary(model)$r.squared
      
    pWilcox[index] = wilcox.test(met ~ urbanRural)$p.value
      
    index=index+1
    
  }
}
  
dFrame <- data.frame(names, meanMetabolite, sdMetabolite, meanUrban, sdUrban, meanRural, sdRural,
                     pValuesUrbanRural, pWilcox)
dFrame$UrbanToRural <- meanUrban/meanRural
dFrame$adjustedPurbanRural <- p.adjust( dFrame$pValuesUrbanRural, method = "BH" )
dFrame$adjustedPwilcox <- p.adjust(dFrame$pWilcox, method="BH")
dFrame$r.squared = r.squared
dFrame <- dFrame[order(dFrame$adjustedPurbanRural),]
write.table(dFrame, file="metabolonModel_pValues.txt", sep="\t",row.names=FALSE, quote=F)

####plots
pathway = read.table("metabolonPathwayInfo.txt", sep="\t", header=T, colClasses="character")
pathway$biochemical = sub("*", "", pathway$biochemical, fixed=T) #* indicates compounds not officially confirmed
##adjust chemical names to have proper punctuation
getMetaboliteName <- function(cheml) {
  if(cheml %in% pathway$biochemical) { #no need to modify
    r = which(pathway$biochemical==cheml)
    return(paste(cheml))
  } else {
    cheml = sub("^X", "", cheml) #need to remove X from chemicals that start with a number
    sp = strsplit(cheml, ".", fixed=T)[[1]] #split on dot, which replaces spaces, dashes, etc
    #find row that matches the split string
    r = grepl(sp[1], pathway$biochemical, fixed=T)
    i=2
    while(length(which(r)) > 1 & i <=length(sp)) {
      r2 = grepl(sp[i], pathway$biochemical, fixed=T)
      r = r & r2
      i = i+1
    }
    if(i > length(sp) & length(which(r)) > 1) { ##match not found
      cheml2 = sub(".", " ", cheml, fixed=T) #androsterone.sulfate has three hits - look for exact match (the other 2 have longer names that match)
      r = cheml2 == pathway$biochemical
      if(length(which(r))==1) {
        return(paste(pathway$biochemical[r]))
      } else if(cheml=="4.androsten.3beta.17beta.diol.disulfate..1..") { #2 hits: [1] "4-androsten-3beta,17beta-diol disulfate (1)" "4-androsten-3beta,17beta-diol disulfate (2)"
        return("4-androsten-3beta,17beta-diol disulfate (1)")
      } else if (cheml=="1.palmitoyl.GPE..16.0.") { #2 hits: [1] "1-(1-enyl-palmitoyl)-GPE (P-16:0)" "1-palmitoyl-GPE (16:0)" 
        return("1-palmitoyl-GPE (16:0)" )
      } else if(cheml=="1.oleoyl.GPC..18.1.") {
        return("1-oleoyl-GPC (18:1)")
      } else if(cheml=="1.stearoyl.GPC..18.0.") {
        return("1-stearoyl-GPC (18:0)")
      } else if(cheml=="pro.hydroxy.pro") {
        return("pro-hydroxy-pro") 
      } else if(cheml=="1.palmitoleoyl.GPC..16.1..") {
        return("1-palmitoleoyl-GPC (16:1)")
      } else if(cheml=="1.palmitoyl.GPC..16.0.") {
        return("1-palmitoyl-GPC (16:0)")
      } else if(cheml=="DSGEGDFXAEGGGVR.") {
        return("DSGEGDFXAEGGGVR")
      } else if(cheml=="1.oleoyl.GPE..18.1.") {
        return("1-oleoyl-GPE (18:1)")
      } else if(cheml=="1.oleoylglycerol..18.1.") {
        return("1-oleoylglycerol (18:1)" )
      } else if(cheml=="1.stearoyl.GPE..18.0.") {
        return("1-stearoyl-GPE (18:0)" )
      } else if(cheml=="beta.alanine") {
        return("beta-alanine")
      } else {
        warning(paste("Multiple hits for chemical:", cheml))
        return(cheml)
      }
    } else {
      return(paste(pathway$biochemical[r]))
    }
  }
}

pdf("metabolon_model_boxplots.pdf")
for(r in 1:nrow(dFrame)) {
  name = dFrame$names[r]
  met = metabolon[,names(metabolon)==name]
  label = paste("metabolite:", getMetaboliteName(name))
  
  graphMain =  paste(label, "\n",
                     " pAdjRuralUrban= ", format(dFrame$adjustedPurbanRural[r],digits=3), sep="")
  boxplot(met ~ urbanRural, main=graphMain, xlab="", ylab="scaled value")
  points(x=urbanRural, y=met, col=ifelse(metabolon$classification == "rural", "blue", "red"), pch=16)
}
dev.off()
