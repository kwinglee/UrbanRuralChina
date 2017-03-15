##generate supplemental tables of significant associations between 16S abundances, metabolite data and metadata.

rm(list=ls())

setwd("")

taxaLevels <- c("phylum","class","order","family","genus", "otu")
levelNames = c("Phylum", "Class", "Order", "Family", "Genus", "OTU")
alpha = 0.05 #adjusted sig level

taxonomy = read.table("abundantOTU.chinaForward.taxonomy.txt", header=F, colClasses="character", sep="\t")
##function that returns full taxonomy for the given name in the given level
getTaxonomy <- function(name, level) {
  if(level == "phylum") {
    return(paste("d__Bacteria;p__", name, sep=""))
  } else if(level=="class") {
    rows = which(taxonomy$V9==name) 
    r = rows[1]
    if(length(unique(taxonomy$V6[rows])) != 1) {
      warning(paste("Mixed class:", name))
    }
    return(paste("d__", taxonomy$V3[r], ";p__", taxonomy$V6[r], ";c__", name, sep=""))
  } else if(level=="order") {
    rows = which(taxonomy$V12==name)
    r = rows[1]
    if(length(unique(taxonomy$V9[rows])) != 1) {
      warning(paste("Mixed order:", name))
    }
    return(paste("d__", taxonomy$V3[r], ";p__", taxonomy$V6[r], ";c__", taxonomy$V9[r], ";o__", name, sep=""))
  } else if(level=="family") {
    rows = which(taxonomy$V15==name)
    r = rows[1]
    if(length(unique(taxonomy$V12[rows])) != 1) {
      ##missing families
      if(name == "Leuconostocaceae") {
        return("d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Leuconostocaceae")
      } else if(name == "Corynebacteriaceae") {
        return("d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Corynebacteriaceae")
      } else if(name == "Clostridiaceae.1") {
        return(getTaxonomy("Clostridiaceae 1", "family"))
      } else if(name == "Oxalobacteraceae") {
        return("d__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Oxalobacteraceae")
      } else {
        warning(paste("Mixed/missing family:", name))
      }
    }
    return(paste("d__", taxonomy$V3[r], ";p__", taxonomy$V6[r], ";c__", taxonomy$V9[r], ";o__", taxonomy$V12[r], ";f__", name, sep=""))
  } else if(level=="genus") {
    if(name=="Escherichia.Shigella") {
      name = "Escherichia/Shigella"
    }
    rows = which(taxonomy$V18==name)
    r = rows[1]
    if(length(unique(taxonomy$V15[rows])) != 1) {
      ##missing genera
      if(name == "Lactococcus") {
        return("d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus")
      } else if(name == "Weissella") {
        return("d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Leuconostocaceae;g__Weissella")
      } else if(name == "Corynebacterium") {
        return("d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Corynebacteriaceae;g__Corynebacterium")
      } else if(name == "Leuconostoc") {
        return("d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Leuconostocaceae;g__Leuconostoc")
      } else if(name == "Atopobium") {
        return("d__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Coriobacteriales;f__Coriobacteriaceae;g__Atopobium")
      } else if(name == "Parvimonas") { 
        return("d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiales Incertae Sedis XI;g__Parvimonas")
      } else if(name == "Cronobacter") {
        return("d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Cronobacter")
      } else {
        warning(paste("Mixed/missing genus:", name))
      }
    }
    return(paste("d__", taxonomy$V3[r], ";p__", taxonomy$V6[r], ";c__", taxonomy$V9[r], ";o__", taxonomy$V12[r], ";f__", taxonomy$V15[r], ";g__", name, sep=""))
  } else if(level=="otu") {
    n = sub("X", "", name)
    n = paste("Consensus", n, sep="")
    r = taxonomy$V1==n
    n = sub("Consensus", "OTU", n)
    return(paste("d__", taxonomy$V3[r], ";p__", taxonomy$V6[r], ";c__", taxonomy$V9[r], ";o__", taxonomy$V12[r], ";f__", taxonomy$V15[r], ";g__", taxonomy$V18[r], ";", n, sep=""))
  } else {
    warning("incorrect level")
  }
}

pathway = read.table("metabolonPathwayInfo.txt", sep="\t", header=T, colClasses="character")
pathway$biochemical = sub("*", "", pathway$biochemical, fixed=T) #* indicates compounds not officially confirmed
##adjust chemical names to have proper punctuation and pathway info
getPathway <- function(cheml) {
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

##function to get metadata labels from given variable var
md.convert = read.table("2015-11-23_metadata//metadataConversion.txt", header=T, sep="\t", colClasses="character")
getMetadataName <- function(var) {
  name = md.convert$Label[md.convert$Variable==var]
  if(is.na(name)) {
    warning(paste("Variable not found:", var))
    return(var)
  } else {
    return(name)
  }
}

##for the given otu table, replace otu name with full taxonomy and add column for level analyzed
formatOTU <- function(taxa, table) {
  lev = rep(taxa, nrow(table))
  name = rep(NA, nrow(table))
  for(i in 1:length(name)) {
    name[i] = getTaxonomy(table$namesOTU[i], taxa)
  }
  return(data.frame(taxLevel=lev, fullName=name))
}

##for the given metadata, replace name with full name
formatMetadata <- function(names) {
  convert = rep(NA, length(names))
  for(i in 1:length(convert)) {
    convert[i] = getMetadataName(names[i])
  }
  return(convert)
}

##for the given metabolon data, replace name with full name
formatMetabolite <- function(names) {
  convert = rep(NA, length(names))
  for(i in 1:length(convert)) {
    convert[i] = getPathway(names[i])
  }
  return(convert)
}

###########metabolite vs microbiome
supptable = data.frame()
for(t in 1:length(taxaLevels)) {
  taxa = taxaLevels[t]
  print(paste("Metabolon vs", taxa))
  ptable = read.table(paste("model_metabolon_v_", taxa, "_pValues.txt", sep=""), header=T, sep="\t")
  sig = ptable[ptable$adjLMp<alpha,]
  print(nrow(sig))
  if(nrow(sig) > 0) {
    sig = sig[order(sig$pLM),]
    format = formatOTU(taxa, sig)
    met = formatMetabolite(sig$namesOther)
    df = data.frame(taxaLevel=rep(levelNames[t], nrow(sig)), 
                    fullName=format$fullName, 
                    pathway=paste(sapply(met, getSuperPathway), 
                                  sapply(met, getSubPathway), 
                                  sep=";"),
                    metabolite=met, 
                    pValue=sig$pLM, 
                    adjustedP=sig$adjLMp,
                    spearman = sig$r,
                    r.squared = sig$r.squared,
                    nonparametricTest = sig$nonparamTest,
                    pNonparametric = sig$pNonparam,
                    pAdjNonparametric = sig$adjNonparmP)
    supptable = rbind(supptable, df)
  }
}
write.table(supptable, "assocTable_metabolite-v-16S.txt", row.names=F, col.names=T, quote=F, sep="\t")


###########metadata vs microbiome
supptable = data.frame() 
for(t in 1:length(taxaLevels)) {
  taxa = taxaLevels[t]
  print(paste("Metadata vs", taxa))
  ptable = read.table(paste("model_metadata_v_", taxa, "_pValues.txt", sep=""), header=T, sep="\t")
  sig = ptable[ptable$adjLMp<alpha,]
  print(nrow(sig))
  if(nrow(sig) > 0) {
    sig = sig[order(sig$pLM),]
    format = formatOTU(taxa, sig)
    df = data.frame(taxaLevel=rep(levelNames[t], nrow(sig)), 
                    fullName=format$fullName, 
                    metadata=formatMetadata(sig$namesOther), 
                    pValue=sig$pLM, 
                    adjustedP=sig$adjLMp,
                    spearman = sig$r,
                    r.squared = sig$r.squared,
                    nonparametricTest = sig$nonparamTest,
                    pNonparametric = sig$pNonparam,
                    pAdjNonparametric = sig$adjNonparmP)
    supptable = rbind(supptable, df)
  }
}
write.table(supptable, "assocTable_metadata-v-16S.txt", row.names=F, col.names=T, quote=F, sep="\t")


###########metadata vs metabolite
supptable = data.frame() 
ptable = read.table(paste("model_metadata_v_metabolon_pValues.txt", sep=""), header=T, sep="\t")
sig = ptable[ptable$adjLMp<alpha,]
print(nrow(sig))
if(nrow(sig) > 0) {
  sig = sig[order(sig$pLM),]
   met = formatMetabolite(sig$namesOTU)
  df = data.frame(pathway=paste(sapply(met, getSuperPathway), 
                                sapply(met, getSubPathway), 
                                sep=";"),
                  metabolite=met, 
                  metadata=formatMetadata(sig$namesOther), 
                  pValue=sig$pLM, 
                  adjustedP=sig$adjLMp,
                  spearman = sig$r,
                  r.squared = sig$r.squared,
                  nonparametricTest = sig$nonparamTest,
                  pNonparametric = sig$pNonparam,
                  pAdjNonparametric = sig$adjNonparmP)
  supptable = rbind(supptable, df)
}
write.table(supptable, "assocTable_metadata-v-metabolite.txt", row.names=F, col.names=T, quote=F, sep="\t")

