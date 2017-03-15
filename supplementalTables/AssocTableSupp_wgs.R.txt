##make supplemental table of significant whole genome sequencing (minikraken and humann)
##associations with metabolite and metadata

rm(list=ls())
setwd("")
meta = c("metabolon", "metadata")
alpha = 0.05 #adjusted sig level

##function to get metadata labels from given variable var
md.convert = read.table("metadataConversion.txt", 
                        header=T, sep="\t", colClasses="character")
getMetadataName <- function(var) {
  name = md.convert$Label[md.convert$Variable==var]
  if(is.na(name)) {
    warning(paste("Variable not found:", var))
    return(var)
  } else {
    return(name)
  }
}

pathway = read.table("metabolonPathwayInfo.txt", 
                     sep="\t", header=T, colClasses="character")
pathway$biochemical = sub("*", "", pathway$biochemical, fixed=T) #* indicates compounds not officially confirmed
##adjust chemical names to have proper punctuation and pathway info
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



###########minikraken
taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")
levelNames = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

##some names need to be cleaned up
clean <- function(name, taxa) {
  if(name == "Deinococcus.Thermus") { #phylum (rest are species)
    return("Deinococcus-Thermus")
  } else if(grepl("X.", name, fixed=T)) { #many of the species start with this
    sp = strsplit(name, "_")[[1]]
    tax = taxa[grepl(paste("_", sp[length(sp)], sep=""), taxa)] #need the underscore or saccharolyticum will match twice
    if(length(tax) != 1) {
      print(tax)
    }
    return(tax)
  } else if(name == "Spiribacter_sp._UAH.SP71") {
    return("Spiribacter_sp._UAH-SP71")
  } else if(name == "Ruminococcus_sp._SR1.5") {
    return("Ruminococcus_sp._SR1/5")
  } else if(grepl("Streptococcus_sp._I", name)) {
    return(sub("_I.", "_I-", name, fixed=T))
  } else if(name == "Vibrio_phage_pYD38.A") {
    return("Vibrio_phage_pYD38-A") 
  } else if(name == "Coprococcus_sp._ART55.1") {
    return("Coprococcus_sp._ART55/1")
  } else if(name == "Anaeromyxobacter_sp._Fw109.5") {
    return("Anaeromyxobacter_sp._Fw109-5")
  } else if(name == "Aster_yellows_witches..broom_phytoplasma") {
    return("Aster_yellows_witches'-broom_phytoplasma")
  } else if(name == "Sulfurovum_sp._NBC37.1") {
    return("Sulfurovum_sp._NBC37-1")
  } else if(name == "Lacinutrix_sp._5H.3.7.4") {
    return("Lacinutrix_sp._5H-3-7-4")
  } else if(name == "Delftia_sp._Cs1.4") {
    return("Delftia_sp._Cs1-4")
  } else if(name == "butyrate.producing_bacterium_SS3.4") {
    return("butyrate-producing_bacterium_SS3/4")
  } else if(name == "butyrate.producing_bacterium_SSC.2") {
    return("butyrate-producing_bacterium_SSC/2")
  } else if(name == "Psychrobacter_sp._PRwf.1") {
    return("Psychrobacter_sp._PRwf-1")
  } else if(name == "Candidatus_Arthromitus_sp._SFB.mouse") {
    return("Candidatus_Arthromitus_sp._SFB-mouse")
  } else if(name == "Pseudovibrio_sp._FO.BEG1") {
    return("Pseudovibrio_sp._FO-BEG1")
  } else if(name == "Agrobacterium_sp._H13.3") {
    return("Agrobacterium_sp._H13-3")
  } else if(name == "Carnobacterium_sp._17.4") {
    return("Carnobacterium_sp._17-4")
  } else if(name == "Oscillatoria_nigro.viridis") {
    return("Oscillatoria_nigro-viridis")
  } else if(name == "Actinoplanes_sp._SE50.110") {
    return("Actinoplanes_sp._SE50/110")
  } else if(name == "Candidatus_Arthromitus_sp._SFB.rat.Yit") {
    return("Candidatus_Arthromitus_sp._SFB-rat-Yit")
  } else if(name == "Candidatus_Sulfuricurvum_sp._RIFRC.1") {
    return("Candidatus_Sulfuricurvum_sp._RIFRC-1")
  } else if(name == "Streptomyces_sp._SirexAA.E") {
    return("Streptomyces_sp._SirexAA-E")
  } else if(name == "Methylobacterium_sp._4.46") {
    return("Methylobacterium_sp._4-46")
  } else if(name == "Spirochaeta_sp._L21.RPul.D2") {
    return("Spirochaeta_sp._L21-RPul-D2")
  } else if(name == "Paenibacillus_sp._JDR.2") {
    return("Paenibacillus_sp._JDR-2")
  } else if(name == "Sphingomonas_sp._MM.1") {
    return("Sphingomonas_sp._MM-1")
  } else if(name == "Dokdonia_sp._4H.3.7.5") {
    return("Dokdonia_sp._4H-3-7-5")
  } else if(name == "Glaciecola_sp._4H.3.7.YE.5") {
    return("Glaciecola_sp._4H-3-7+YE-5")
  } else if(name == "Enterobacteria_phage_fiAA91.ss") {
    return("Enterobacteria_phage_fiAA91-ss")
  } else if(name == "Bacteroides_phage_B124.14") {
    return("Bacteroides_phage_B124-14")
  } else if(name == "Pandoraea_sp._RB.44") {
    return("Pandoraea_sp._RB-44")
  } else if(name == "Actinoplanes_sp._N902.109") {
    return("Actinoplanes_sp._N902-109")
  } else if(name == "Streptococcus_phage_YMC.2011") {
    return("Streptococcus_phage_YMC-2011")
  } else if(name == "Sphingobium_sp._SYK.6") {
    return("Sphingobium_sp._SYK-6")
  } else if(name == "Pantoea_sp._At.9b") {
    return("Pantoea_sp._At-9b")
  } else if(name == "Roseiflexus_sp._RS.1") {
    return("Roseiflexus_sp._RS-1")
  } else if(name == "Synechococcus_sp._JA.3.3Ab") {
    return("Synechococcus_sp._JA-3-3Ab")
  } else if(name == "Synechococcus_sp._JA.2.3B.a.2.13.") {
    return("Synechococcus_sp._JA-2-3B'a(2-13)")
  } else if(name == "Yersinia_phage_L.413C") {
    return("Yersinia_phage_L-413C")
  } else if(name == "Pusillimonas_sp._T7.7") {
    return("Pusillimonas_sp._T7-7")
  } else if(name == "Shewanella_sp._W3.18.1") {
    return("Shewanella_sp._W3-18-1")
  } else if(name == "Flavobacteriaceae_bacterium_3519.10") {
    return("Flavobacteriaceae_bacterium_3519-10")
  } else if(name == "Enterobacter_sp._R4.368") {
    return("Enterobacter_sp._R4-368")
  } else {
    return(name)
  }
}

###minikraken vs metadata/metabolon
for(m in meta) {
  print(m)
  supptable = data.frame()
  for(i in 1:length(taxaLevels)) {
    lev = taxaLevels[i]
    print(lev)
    table = read.table(paste("minikraken_model_pValues_", m, "_v_", lev, ".txt", sep=""),
                       header=T, sep='\t')
    taxonomy = read.table(paste("minikraken_merged_", lev, ".txt", sep=""),
                          header=T, sep="\t", stringsAsFactors = F, quote="")
    ##subset to significant columns only
    sig = table[table$adjLMp<alpha,]
    print(nrow(sig))
    if(nrow(sig) > 0) {
      sig = sig[order(sig$pLM),]
      
      ##get full taxonomy and other name
      fulltaxa = rep(NA, nrow(sig))
      other = rep(NA, nrow(sig))
      metPath = rep(NA, nrow(sig))
      for(r in 1:nrow(sig)) {
        name = clean(sig$namesOTU[r], taxonomy$taxa)
        fulltaxa[r] = taxonomy$taxonomy[taxonomy$taxa==name]
        if(m == "metabolon") {
          other[r] = getMetaboliteName(table$namesOther[r])
          metPath[r] = paste(getSuperPathway(other[r]), 
                             getSubPathway(other[r]), 
                             sep=";")
        } else {
          other[r] = getMetadataName(table$namesOther[r])
        }
      }
      
      if(m == "metabolon") {
        df = data.frame(taxaLevel=rep(levelNames[i], nrow(sig)), 
                        fullName=fulltaxa, 
                        metPath,
                        other, 
                        pLM=sig$pLM, 
                        pAdjLM=sig$adjLMp,
                        spearman = sig$r,
                        r.sqaured = sig$r.squared,
                        nonparametricTest = sig$nonparamTest,
                        pNonparametric = sig$pNonparam,
                        pAdjNonparametric = sig$adjNonparmP)
        
      } else {
        df = data.frame(taxaLevel=rep(levelNames[i], nrow(sig)), 
                        fullName=fulltaxa, 
                        other, 
                        pLM=sig$pLM, 
                        pAdjLM=sig$adjLMp,
                        spearman = sig$r,
                        r.sqaured = sig$r.squared,
                        nonparametricTest = sig$nonparamTest,
                        pNonparametric = sig$pNonparam,
                        pAdjNonparametric = sig$adjNonparmP)
      }
      
      supptable = rbind(supptable, df)
    }
  }
  write.table(supptable, paste("assocTable_", m, "-v-minikraken.txt", sep=""), 
              row.names=F, col.names=T, quote=F, sep="\t")
}


###########humann
levels = c("module", "pathway")
levelNames = c("Module", "Pathway")

###humann vs metadata/metabolon
for(m in meta) {
  print(m)
  supptable = data.frame()
  for(i in 1:length(levels)) {
    lev = levels[i]
    print(lev)
    table = read.table(paste("humann_model_unlog_", m, "_v_", lev, "_pValues.txt", sep=""),
                       header=T, sep='\t', comment.char="", quote="")
    ##subset to significant columns only
    sig = table[table$adjLMp<alpha,]
    print(nrow(sig))
    if(nrow(sig) > 0) {
      sig = sig[order(sig$pLM),]
      
      ##correct names
      other = rep(NA, nrow(sig))
      if(m == "metabolon") {
        for(r in 1:nrow(sig)) {
          other[r] = getMetaboliteName(table$namesOther[r])
        } 
        df = data.frame(keggLevel=rep(levelNames[i], nrow(sig)), 
                        fullName=sig$description, 
                        pathway=paste(sapply(other, getSuperPathway), 
                                      sapply(other, getSubPathway), 
                                      sep=";"),
                        other, 
                        pLM=sig$pLM, 
                        pAdjLM=sig$adjLMp,
                        spearman = sig$r,
                        r.squared = sig$r.squared,
                        nonparametricTest = sig$nonparamTest,
                        pNonparametric = sig$pNonparam,
                        pAdjNonparametric = sig$adjNonparmP)
      } else {
        for(r in 1:nrow(sig)) {
          other[r] = getMetadataName(table$namesOther[r])
        }
        
        df = data.frame(keggLevel=rep(levelNames[i], nrow(sig)), 
                        fullName=sig$description, 
                        other, 
                        pLM=sig$pLM, 
                        pAdjLM=sig$adjLMp,
                        spearman = sig$r,
                        r.squared = sig$r.squared,
                        nonparametricTest = sig$nonparamTest,
                        pNonparametric = sig$pNonparam,
                        pAdjNonparametric = sig$adjNonparmP)
      }
      supptable = rbind(supptable, df)
    }
  }
  write.table(supptable, paste("assocTable_", m, "-v-humann.txt", sep=""), 
              row.names=F, col.names=T, quote=F, sep="\t")
}