##Plot prevalence vs. p-value on rolling window

rm(list=ls())
setwd("")

table = read.table("rollingWindownPrevelance.txt", header = T, sep='\t', 
                   colClasses = c(rep("numeric", 3), "logical"))

tiff("rollingPrevalenceVp.tif", res=300, height=900, width=1600)
par(mar=c(4.1, 4.1, .5, 7.1))
plot(x=table$pValue[table$higherInRural], y=table$rollingPrevelance[table$higherInRural],
     xlab="unadjusted P-value", ylab="prevalence", log="x",
     type="l", col="blue", ylim=range(table$rollingPrevelance))
lines(x=table$pValue[!table$higherInRural], y=table$rollingPrevelance[!table$higherInRural],
      col="red")
legend("topright", xpd=T, inset=c(-.45, 0),
       legend=c("higher in rural", "higher in urban"),
       col=c("blue", "red"), lty=c(1,1), cex=.8)
dev.off()