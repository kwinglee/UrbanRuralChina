##compare non-parametric and linear model results

rm(list=ls())
setwd("")

##return number significant in linear model, non parametric, and both for individual taxa/genes/metadata/metabolites
otuModel <- function(dataset, level=NA) {
  table = read.table(paste("SuppTableModel_", dataset, ifelse(is.na(level), "", paste("_", level, sep="")),
                           ".txt", sep=""), stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")
  if("pAdjUrbanRural" %in% names(table)) {
    names(table)[names(table)=="pAdjUrbanRural"] = "pAdjUrbanRuralLM"
  }
  if("pAdjUrbanRuralNonparametric" %in% names(table)) {
    names(table)[names(table)=="pAdjUrbanRuralNonparametric"] = "pAdjUrbanRuralWilcox"
  }
  return(c(sum(table$pAdjUrbanRuralLM < 0.05, na.rm=T), sum(table$pAdjUrbanRuralWilcox < 0.05),
           sum(table$pAdjUrbanRuralLM < 0.05 & table$pAdjUrbanRuralWilcox < 0.05)))
}

##PCoA model return if MDS1 significant in linear model, non parametric, and both
pcoaModel <-function(tableName) {
  table=read.table(tableName, stringsAsFactors = F, sep="\t", header = T)
  if("pValuesRuralUrban" %in% names(table)) {
    names(table)[names(table=="pValuesRuralUrban")] = "pValuesUrbanRural"
  }
  if("pValuesRuralUrbanWilcox" %in% names(table)) {
    names(table)[names(table=="pValuesRuralUrbanWilcox")] = "pValuesUrbanRuralWilcox"
  }
  return(c(table$pValuesUrbanRural[table$names=="MDS1"] < 0.05,
         table$pValuesUrbanRuralWilcox[table$names=="MDS1"] < 0.05,
         table$pValuesUrbanRural[table$names=="MDS1"] < 0.05 & table$pValuesUrbanRuralWilcox[table$names=="MDS1"] < 0.05))
}

names = c("data", "level", "numLinearModel", "numNonParam", "numBoth", "MDS1linear", "MDS1nonParam", "MDS1both")
results = data.frame(data=character(), level=character(), numLinearModel=numeric(),
                     numNonParam=numeric(), numBoth=numeric(), MDS1linear=vector(),
                     MDS1nonParam=vector(), MDS1both=vector(), stringsAsFactors = F)
row=1

##metadata
otu = otuModel("metadata")
pc = pcoaModel("pcoaModel_pValues_metadata.txt")
df = c("metadata", NA, otu, pc)
results[row,] = df
row=row+1

##metabolite
otu = otuModel("metabolome")
pc = pcoaModel("pcoaModel_pValues_metabolon.txt")
df = c("metabolome", NA, otu, pc)
results[row,] = df
row=row+1

##humann
levels = c("module", "pathway")
for(lev in levels) {
  otu = otuModel("humann", lev)
  pc = pcoaModel(paste("humann_pcoaModel_pValues_", lev, ".txt", sep=""))
  results[row,] = c("gene WGS", lev, otu, pc)
  row=row+1
}

##minikraken
taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")
for(lev in taxaLevels) {
  otu = otuModel("minikraken", lev)
  pc = pcoaModel(paste("minikraken_pcoaModel_pValues_", lev, ".txt", sep=""))
  results[row,] = c("taxonomic WGS", lev, otu, pc)
  row=row+1
}


names(results) = names
write.table(results, "lm_vs_nonparametric.txt", row.names = F, col.names = T, quote = F, sep="\t")