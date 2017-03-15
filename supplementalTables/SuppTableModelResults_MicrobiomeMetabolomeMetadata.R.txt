##supplemental table of model results for 16S, metabolite and metadata

rm(list=ls())
setwd("")

####microbiome
taxaLevels = c("phylum", "class", "order", "family", "genus", "otu")

taxonomy = read.table("RDPtaxonomy.txt", header=T, colClasses="character", sep="\t")
otu.tax = read.table("abundantOTU.chinaForward.taxonomy.txt", header=F, colClasses="character", sep="\t")
##function that returns full taxonomy for the given name in the given level
getTaxonomy <- function(name, level) {
  name = gsub(".", " ", name, fixed=T)
  if(level == "phylum") {
    return(paste("d__Bacteria;p__", name, sep=""))
  } else if(level=="class") {
    rows = which(taxonomy$class==name) 
    r = rows[1]
    if(length(unique(taxonomy$class[rows])) != 1) {
      warning(paste("Mixed class:", name))
    }
    return(paste("d__Bacteria;p__", taxonomy$phylum[r], ";c__", name, sep=""))
  } else if(level=="order") {
    rows = which(taxonomy$order==name)
    r = rows[1]
    if(length(unique(taxonomy$order[rows])) != 1) {
      warning(paste("Mixed order:", name))
    }
    return(paste("d__Bacteria;p__", taxonomy$phylum[r], ";c__", taxonomy$class[r], ";o__", name, sep=""))
  } else if(level=="family") {
    rows = which(taxonomy$family==name)
    r = rows[1]
    if(length(unique(taxonomy$family[rows])) != 1) {
      warning(paste("Mixed/missing family:", name))
    }
    return(paste("d__Bacteria;p__", taxonomy$phylum[r], ";c__", taxonomy$class[r], ";o__", taxonomy$order[r], ";f__", name, sep=""))
  } else if(level=="genus") {
    rows = which(taxonomy$genus==name)
    r = rows[1]
    if(length(unique(taxonomy$genus[rows])) != 1) {
      if(name == "Escherichia Shigella") {
        return(getTaxonomy("Escherichia/Shigella", "genus"))
      } else if(name == "Weissella") {
        return("d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Leuconostocaceae;g__Weissella")
      } else {
        warning(paste("Mixed/missing genus:", name))
      }
    } else {
      return(paste("d__Bacteria;p__", taxonomy$phylum[r], ";c__", taxonomy$class[r], ";o__", taxonomy$order[r], ";f__", taxonomy$family[r], ";g__", name, sep=""))  
    }
  } else if(level=="otu") {
    n = sub("OTU", "Consensus", name)
    r = otu.tax$V1==n
    return(paste("d__Bacteria;p__", otu.tax$V6[r], ";c__", otu.tax$V9[r], ";o__", otu.tax$V12[r], ";f__", otu.tax$V15[r], ";g__", otu.tax$V18[r], ";", name, sep=""))
  } else {
    warning("incorrect level")
  }
}

for(taxa in taxaLevels) {
  table = read.table(paste("otuModel_pValues_", taxa, ".txt", sep=""), header=T, sep="\t", colClasses=c("character", rep("numeric", 17)))
  
  if(taxa == "otu") {
    table$names = sub("X", "OTU", table$names)
  }
  
  df = data.frame(taxa = table$names, 
                  fullTaxonomy = sapply(table$names, getTaxonomy, level=taxa, USE.NAMES = F),
                  meanRuralT1 = table$meanRural1,
                  sdRuralT1 = table$sdRural1,
                  meanRuralT2 = table$meanRural2,
                  sdRuralT2 = table$sdRural2,
                  meanRural = table$meanRural,#(table$meanRural1+table$meanRural2)/2,
                  sdRural = table$sdRural,
                  meanUrban1 = table$meanUrban1,
                  sdUrban1 = table$sdUrban1,
                  meanUrban2 = table$meanUrban2,
                  sdUrban2 = table$sdUrban2,
                  meanUrban = table$meanUrban,#(table$meanUrban1+table$meanUrban2)/2,
                  sdUrban = table$sdUrban,
                  pUrbanRural = table$pValuesUrbanRural,
                  pTime = table$pValuesTime,
                  pSubj = table$pValuesSubject,
                  pAdjUrbanRural = table$adjustedPurbanRural,
                  pAdjTime = table$adjustedPtime,
                  pAdjSubj = table$adjustedPsubject,
                  rSquared = table$r.squared,
                  pUrbanRuralWilcoxT1 = table$pValuesUrbanRuralWilcoxT1,
                  pUrbanRuralWilcoxT2 = table$pValuesUrbanRuralWilcoxT2,
                  pUrbanRuralWilcoxAllTimes = table$pValuesUrbanRuralWilcoxAll,
                  pAdjUrbanRuralWilcoxT1 = table$adjustedPurbanRuralWilcoxT1,
                  pAdjUrbanRuralWilcoxT2 = table$adjustedPurbanRuralWilcoxT2,
                  pAdjUrbanRuralWilcoxAllTimes = table$adjustedPurbanRuralWilcoxAll)
  df = df[order(df$pUrbanRural),]
  
  write.table(df, paste("SuppTableModel_16S_", taxa, ".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
}

####metabolome
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

##return the super pathway for the given cleaned up chemical name
getSuperPathway <- function(cheml) {
  r = which(pathway$biochemical==cheml)
  if(length(r) != 1) {
    warning(paste("Multiple hits for chemical:", cheml))
  }
  return(pathway$super.pathway[r])
}

##return the sub pathway for the given cleaned up chemical name
getSubPathway <- function(cheml) {
  r = which(pathway$biochemical==cheml)
  if(length(r) != 1) {
    warning(paste("Multiple hits for chemical:", cheml))
  }
  return(pathway$sub.pathway[r])
}

table = read.table("metabolonModel_pValues.txt", header=T, sep="\t", colClasses=c("character", rep("numeric", 9)))
met = sapply(table$names, getMetaboliteName, USE.NAMES=F)
df = data.frame(metabolite = met,
                superPathway = sapply(met, getSuperPathway, USE.NAMES=F),
                subPathway = sapply(met, getSubPathway, USE.NAMES=F),
                meanRural = table$meanRural,
                sdRural = table$sdRural,
                meanUrban = table$meanUrban,
                sdUrban = table$sdUrban,
                pUrbanRuralLM = table$pValuesUrbanRural,
                pAdjUrbanRuralLM = table$adjustedPurbanRural,
                rSquared = table$r.squared,
                pUrbanRuralWilcox = table$pWilcox,
                pAdjUrbanRuralWilcox = table$adjustedPwilcox)
df = df[order(df$pUrbanRuralLM),]
write.table(df, "SuppTableModel_metabolome.txt", sep="\t", quote=F, row.names=F, col.names=T)

####metadata
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

table = read.table("metadataModel_pValues.txt", header=T, sep="\t")
df = data.frame(metadata = sapply(table$names, getMetadataName, USE.NAMES=F),
                meanRural = table$meanRural,
                sdRural = table$sdRural,
                meanUrban = table$meanUrban,
                sdUrban = table$sdUrban,
                pUrbanRuralLM = table$pValuesUrbanRural,
                pAdjUrbanRural = table$adjustedPurbanRural,
                rSquared = table$r.squared,
                nonparametricTest = table$nonparametricTest,
                pUrbanRuralNonparametric = table$pNonparametric,
                pAdjUrbanRuralNonparametric = table$adjustedPnonparametric)
df = df[order(df$pUrbanRuralLM),]
write.table(df, "SuppTableModel_metadata.txt", sep="\t", quote=F, row.names=F, col.names=T)
