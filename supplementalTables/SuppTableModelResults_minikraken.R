##make supplemental table of whole genome sequencing minikraken model results

rm(list=ls())
setwd("")

taxaLevels = c("domain", "phylum", "class", "order", "family", "genus", "species")

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

for(lev in taxaLevels) {
  print(lev)
  table = read.table(paste("minikraken_otuModel_pValues_", lev, ".txt", sep=""),
                     header=T, sep="\t", stringsAsFactors = F)
  taxonomy = read.table(paste("minikraken_merged_", lev, ".txt", sep=""),
                        header=T, sep="\t", stringsAsFactors = F, quote="")
  fulltaxa = rep(NA, nrow(table))
  for(i in 1:nrow(table)) {
    table$names[i] = clean(table$names[i], taxonomy$taxa)
    fulltaxa[i] = taxonomy$taxonomy[taxonomy$taxa==table$names[i]]
  }
  if(any(duplicated(table$names))) {
    stop(paste("Duplicated name:", table$names[duplicated(table$names)]))
  }
  
  df = data.frame(taxa = table$names,
                  fullTaxonomy = fulltaxa,
                  meanRural = table$meanRural,
                  sdRural = table$sdRural,
                  meanUrban = table$meanUrban,
                  sdUrban = table$sdUrban,
                  pUrbanRuralLM = table$pValuesUrbanRural,
                  pAdjUrbanRuralLM = table$adjustedPurbanRural,
                  rSquared = table$r.squared,
                  pUrbanRuralWilcox = table$pValuesUrbanRuralWilcox,
                  pAdjUrbanRuralWilcox = table$adjustedPurbanRuralWilcox)
  df = df[order(df$pUrbanRuralLM),]
  write.table(df, paste("SuppTableModel_minikraken_", lev, ".txt", sep=""), sep="\t",
              row.names = F, col.names = T, quote=F)
}