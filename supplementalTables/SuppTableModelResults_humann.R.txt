##make supplemental table of whole genome sequencing humann model results

rm(list=ls())

setwd("")

levels = c("module", "pathway")

for(lev in levels) {
  print(lev)
  table = read.table(paste("humann_otuModel_pValues_unlog_", lev, ".txt",sep=""), sep="\t", header=T)
  
  df = data.frame(kegg = table$kegg,
                  description = table$name,
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
  write.table(df, paste("SuppTableModel_humann_", lev, ".txt", sep=""), sep="\t",
              row.names = F, col.names = T, quote=F)
}
