##make tables of significant associations with proportion of reads aligned to MvirDB (virulence genes)

rm(list=ls())
setwd("")
alpha = 0.05 #adjusted sig level

#####16S
taxaLevels <- c("phylum","class","order","family","genus")
levelNames = c("Phylum", "Class", "Order", "Family", "Genus")

taxonomy = read.table("abundantOTU.chinaForward.taxonomy.txt", header=F, colClasses="character", sep="\t")
##function that returns full taxonomy for the given name in the given level
get16sTaxonomy <- function(name, level) {
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

supptable = data.frame()
for(i in 1:length(taxaLevels)) {
  lev = taxaLevels[i]
  print(lev)
  table = read.table(paste("virulenceAssoc_pValues_MvirDB_v_16S_", lev, ".txt", sep=""),
                     header=T, sep="\t")
  sig = table[table$adjLMp < alpha,]
  print(nrow(sig))
  if(nrow(sig) > 0) {
    sig = sig[order(sig$pLM),]
    fulltax = sapply(as.character(sig$namesOther), get16sTaxonomy, lev, USE.NAMES = F)
    df = data.frame(taxaLevel=rep(levelNames[i], nrow(sig)), 
                    fullName=fulltax, 
                    pLM=sig$pLM, 
                    pAdjLM=sig$adjLMp,
                    spearman=sig$r,
                    r.squared=sig$r.squared,
                    pKendall=sig$kendallP,
                    pAdjKendall=sig$adjKendallP)
    supptable = rbind(supptable, df)
  }
}
##otu different format
table = read.table(paste("virulenceAssoc_pValues_MvirDB_v_16S_otu.txt", sep=""),
                   header=T, sep="\t")
sig = table[table$adjLMp < alpha,]
print(nrow(sig))
sig = sig[order(sig$pLM),]
sig$description = gsub("_", "__", as.character(sig$description))
sig$description = sub("Consensus", "OTU", sig$description)
df = data.frame(taxaLevel=rep("OTU", nrow(sig)), 
                fullName=paste("d__Bacteria;", sig$description, sep=""),
                pLM=sig$pLM, 
                pAdjLM=sig$adjLMp,
                spearman=sig$r,
                r.squared=sig$r.squared,
                pKendall=sig$kendallP,
                pAdjKendall=sig$adjKendallP)
supptable = rbind(supptable, df)
write.table(supptable, "virAssocTable_16S.txt", sep="\t", row.names = F, col.names = T, quote=F)

#####metabolon
table = read.table("virulenceAssoc_pValues_MvirDB_v_metabolon.txt",
                   header=T, sep="\t")
sig = table[table$adjLMp < alpha,] #no hits so don't do anything
print(nrow(sig))

#####metadata
table = read.table("virulenceAssoc_pValues_MvirDB_v_metadata.txt",
                   header=T, sep="\t")
sig = table[table$adjLMp < alpha,] #no hits so don't do anything
print(nrow(sig))

#####minikraken
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

supptable = data.frame()
for(i in 1:length(taxaLevels)) {
  lev = taxaLevels[i]
  print(lev)
  table = read.table(paste("virulenceAssoc_pValues_MvirDB_v_minikraken_", lev, ".txt", sep=""),
                     header=T, sep="\t")
  taxonomy = read.table(paste("minikraken_merged_", lev, ".txt", sep=""),
                        header=T, sep="\t", stringsAsFactors = F, quote="")
  sig = table[table$adjLMp < alpha,]
  print(nrow(sig))
  if(nrow(sig) > 0) {
    sig = sig[order(sig$pLM),]
    names = sapply(as.character(sig$namesOther), clean, taxonomy, USE.NAMES = F)
    fulltax = sapply(names, function(x){return(taxonomy$taxonomy[taxonomy$taxa==x])}, USE.NAMES = F)
    df = data.frame(taxaLevel=rep(levelNames[i], nrow(sig)), 
                    fullName=fulltax, 
                    pLM=sig$pLM, 
                    pAdjLM=sig$adjLMp,
                    spearman=sig$r,
                    r.squared=sig$r.squared,
                    pKendall=sig$kendallP,
                    pAdjKendall=sig$adjKendallP)
    supptable = rbind(supptable, df)
  }
}
write.table(supptable, "virAssocTable_minikraken.txt", sep="\t", row.names = F, col.names = T, quote=F)

#####humann
levels = c("module", "pathway")
levelNames = c("Module", "Pathway")
supptable = data.frame()
for(i in 1:length(levels)) {
  lev = levels[i]
  print(lev)
  table = read.table(paste("virulenceAssoc_pValues_MvirDB_v_humann_", lev, ".txt", sep=""),
                     header=T, sep="\t")
  sig = table[table$adjLMp < alpha,]
  print(nrow(sig))
  if(nrow(sig) > 0) {
    sig = sig[order(sig$pLM),]
    df = data.frame(level=rep(levelNames[i], nrow(sig)), 
                    fullName=sig$description, 
                    pLM=sig$pLM, 
                    pAdjLM=sig$adjLMp,
                    spearman=sig$r,
                    r.squared=sig$r.squared,
                    pKendall=sig$kendallP,
                    pAdjKendall=sig$adjKendallP)
    supptable = rbind(supptable, df)
  }
}
write.table(supptable, "virAssocTable_humann.txt", sep="\t", row.names = F, col.names = T, quote=F)
