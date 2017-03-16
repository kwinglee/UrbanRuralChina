##model results for metadata data using linear models and nonparametric (Wilcox) test

rm(list=ls())

setwd("")

metadata = read.table("metadata_cleaned11-2015.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 55)))

categoricalVars = c("F9", "F13", "F17", "sex") #categorical variables

###p-values
##output vectors
names <- vector()
pValuesUrbanRural <- vector()
pNonparametric <- vector()
nonparametricTest <- vector()
mean <- vector()
sd <- vector()
meanUrban <- vector()
sdUrban <- vector()
meanRural <- vector()
sdRural <- vector()
r.squared <- vector()
index <- 1
ruralUrban <- factor(metadata$ruralUrban)

for(i in 3:ncol(metadata)) {
  met <- metadata[,i]
  mean[index] <- mean(met)
  meanUrban[index] <- mean(met[metadata$ruralUrban=="urban"])
  meanRural[index] <- mean(met[metadata$ruralUrban=="rural"])
  sd[index] <- sd(met)
  sdUrban[index] <- sd(met[metadata$ruralUrban=="urban"])
  sdRural[index] <- sd(met[metadata$ruralUrban=="rural"])
  names[index] = names(metadata)[i]
  
  if(names(metadata)[i] %in% categoricalVars) {
    groups = sort(unique(met))
    if(length(groups) != 2) {
      warning(paste("Bad category", i))
    }
    pNonparametric[index] = fisher.test(x=ruralUrban, y=met)$p.value
    nonparametricTest[index] = "Fisher"
    
    pValuesUrbanRural[index] = NA
  } else {
    pNonparametric[index] = wilcox.test(met ~ ruralUrban)$p.value
    nonparametricTest[index] = "Wilcox"
    
    ##linear model
    model = lm(met ~ ruralUrban)
    
    pValuesUrbanRural[index] = anova(model)$`Pr(>F)`[1]
    
    r.squared[index] = summary(model)$r.squared
  }
  
  index=index+1
  
}

dFrame <- data.frame(names, mean, sd, meanUrban, sdUrban, meanRural, sdRural,
                     pValuesUrbanRural, nonparametricTest, pNonparametric)
dFrame$UrbanToRural <- meanUrban/meanRural
dFrame$adjustedPurbanRural <- p.adjust(dFrame$pValuesUrbanRural, method = "BH")
dFrame$adjustedPnonparametric<- p.adjust(dFrame$pNonparametric, method="BH")
dFrame$r.squared = r.squared
dFrame <- dFrame [order(dFrame$adjustedPurbanRural),]
write.table(dFrame, file="metadataModel_pValues.txt", sep="\t",row.names=FALSE, quote=F)

###plots
##function to get metadata labels from given variable var
md.convert = read.table("metadataConversion.txt", header=T, sep="\t", colClasses="character")
getMetadataName <- function(var) {
  name = md.convert$Label[md.convert$Variable==var]
  if(is.na(name)) {
    warning(paste("Variable not found:", var))
    return(var)
  } else {
    return(name)
  }
}

pdf("metadata_model_boxplots.pdf")
for(r in 1:nrow(dFrame)) {
  name = dFrame$names[r]
  met = metadata[,names(metadata)==name]
  label = paste("metadata:", getMetadataName(name))
  
  graphMain =  paste(label, "\n",
                     " pAdjRuralUrban= ", format(dFrame$adjustedPurbanRural[r],digits=3), sep="")
  boxplot(met ~ ruralUrban, main=graphMain, xlab="", ylab=getMetadataName(name))
  points(x=ruralUrban, y=met, col=ifelse(metadata$ruralUrban == "rural", "blue", "red"), pch=16)
}
dev.off()
