##Generate supplemental figure of correlations between variables measured in both metadata and metabolite 

rm(list=ls())

setwd("")

metadata = read.table("metadata_cleaned11-2015.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 55)))
metabolon = read.table("MetabolonScaledImpData-transpose.txt", header=T, sep="\t", colClasses=c("numeric", "character", rep("numeric", 337)))
pValues = read.table("model_metadata_v_metabolon_pValues.txt", header=T, sep="\t", colClasses=c("character", "character", rep("numeric", 5)))

##get sample ids to be the same and merge
names(metabolon)[1] = "sampleID"
names(metabolon)[2] = "ruralUrban"
mg = merge(metadata, metabolon, by=c("sampleID", "ruralUrban"))

##set up plots
col = ifelse(mg$ruralUrban == "rural", "blue", "red")

tiff("sameVarMetadataMetaboliteCorr.tif", res=300, height=4200, width=8000)
par(mfrow=c(1,2), cex.lab=1.5, cex=2.5, cex.main=1.5, mar=c(5, 5, 5, 2), oma=c(2,0,0,0))

##cholesterol
plot(x = mg$cholesterol, y = mg$CHOL, 
     main = paste("Cholesterol\nAdjusted P = ", 
                  format(pValues$adjLMp[pValues$namesOTU=="cholesterol" & pValues$namesOther=="CHOL"], digits=2), 
                  "\nSpearman Correlation = ", 
                  format(pValues$r[pValues$namesOTU=="cholesterol" & pValues$namesOther=="CHOL"], digits=2), sep=""), 
     xlab="metabolite data (scaled abundance)", ylab="metadata (mmol/L)", 
     pch=16, col = col)
abline(lm(mg$CHOL~mg$cholesterol))
mtext("A", side=3, line=2.5, adj=0, cex=5)

##glucose
plot(x = mg$glucose, y = mg$GLU, 
     main = paste("Glucose\nAdjusted P = ", 
                  format(pValues$adjLMp[pValues$namesOTU=="glucose" & pValues$namesOther=="GLU"], digits=2), 
                  "\nSpearman Correlation = ", 
                  format(pValues$r[pValues$namesOTU=="glucose" & pValues$namesOther=="GLU"], digits=2), sep=""), 
     xlab="metabolite data (scaled abundance)", ylab="metadata (mmol/L)", 
     pch=16, col = col)
abline(lm(mg$GLU~mg$glucose))
mtext("B", side=3, line=2.5, adj=0, cex=5)

##legend
par(oma=c(0,0,0,0), mar=c(0.1,0,0,0), new=T, xpd=T, fig=c(0,1,0,1))
plot(1, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")
leg = legend("bottom", horiz=T,
       legend=c("rural", "urban", "linear regression line"),
       col=c("blue", "red", "black"),
       pch=c(16, 16, NA), 
       lty=c(NA, NA, 1),
       cex=1,
       x.intersp = .5,
       text.width=c(.05, .05, .1),
       bty="n")
rect(leg$rect$left[3], leg$rect$top-leg$rect$h, leg$rect$left[1]+leg$rect$w[3], leg$rect$h-.04, border="black")

dev.off()
