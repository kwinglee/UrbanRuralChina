##compare size of number of genes in closest NCBI genome

rm(list=ls())
setwd("")

ptab = read.table("otuModel_pValues_otu.txt", header=T, sep="\t", stringsAsFactors = F)
nctab = read.table("ncbiBlastResultsWithGeneStats.txt", header=T, sep="\t", stringsAsFactors = F)

##subset to just the significant unadjusted and check numbers
ptab = ptab[ptab$pValuesUrbanRural < 0.05,c(1,4,10)] #names, meanUrban and meanRural
ptab$names = sub("X", "Consensus", ptab$names)
nrow(ptab[ptab$meanUrban < ptab$meanRural,]) #80
nrow(ptab[ptab$meanUrban > ptab$meanRural,]) #91
names(ptab)[1] = "OTU"

##merge
mrg = merge(ptab, nctab, by="OTU")

hiRur = mrg[mrg$meanUrban < mrg$meanRural,6]
hiUrb = mrg[mrg$meanUrban > mrg$meanRural,6]
df = rbind(data.frame(group=rep("higher in rural", length(hiRur)), numGenes=hiRur),
           data.frame(group=rep("higher in urban", length(hiUrb)), numGenes=hiUrb))

##plot
p = t.test(df$numGenes~df$group)$p.value
print(p)
tiff("ncbiGeneNumber.tif", res=300, height=1500, width=1200)
par(mar=c(3,4,4,1))
boxplot(df$numGenes~df$group,
        ylab="number of genes",
        main=paste("NCBI Bacterial Genomes\np = ", format(p, digits = 2)),
        pch=16)
dev.off()