##heatmap figure of means of significant OTUs (all taxonomic levels) for rural/urban
##combined with database boxplots

rm(list=ls())
setwd("")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(cowplot)
library(gtable)

taxonomy = read.table("abundantOTU.chinaForward.taxonomy.txt", header=F, colClasses="character", sep="\t")
##function that returns full taxonomy for the given name in the given level
getTaxonomy <- function(name, level) {
  if(level=="otu") {
    num = sub("X", "", name)
    n = paste("Consensus", num, sep="")
    r = taxonomy$V1==n
    ret = paste(taxonomy$V18[r], " (OTU", num, ")", sep="")
    ret = gsub("_", " ", ret) #remove underscores from Lachno
    return(ret)
  } else {
    return(name)
  }
}

####
##combine phylum, family, genus (class and order had nothing significant)
taxaLevels = c("phylum", "family", "genus")

##function to return data frame with means for significant OTUs
getDataFrame <- function(taxa) {
  print(taxa)
  inFileName <- paste("otuModel_pValues_", taxa,  ".txt", sep ="")
  table <-read.table(inFileName,header=TRUE,sep="\t")
  numCols <- ncol(table)
  ColClasses <- c("character", rep("numeric", numCols-1))
  table <-read.table(inFileName,header=TRUE,sep="\t",colClasses=ColClasses)
  
  ##sort by pUrbanRural
  table = table[order(table$pValuesUrbanRural),]
  
  sig = table$adjustedPurbanRural < 0.05
  df = data.frame(names=table$names[sig],meanRural1=table$meanRural1[sig], meanRural2=table$meanRural2[sig],
                  meanUrban1=table$meanUrban1[sig], meanUrban2=table$meanUrban2[sig])
  df$names = as.character(df$names)
  for(i in 1:nrow(df)) {
    df$names[i] = getTaxonomy(df$names[i], taxa)
  } 
  return(rbind(data.frame(names=taxa, meanRural1=NA, meanRural2=NA, meanUrban1=NA, meanUrban2=NA), df))
  
}

df1 = data.frame(names=character(), meanRural1=numeric(), meanRural2=numeric(), meanUrban1=numeric(), meanUrban2=numeric())

for(taxa in taxaLevels) {
  df1 = rbind(df1, getDataFrame(taxa))
}
##give different number of spaces for blank rows so when melt will come out correctly
df1$names = as.character(df1$names)
df1$names[df1$names=="phylum"] = " "
df1$names[df1$names=="family"] = "  "
df1$names[df1$names=="genus"] = "   "

########
##otu
df2 = getDataFrame("otu")
df2$names = as.character(df2$names)
df2$names[df2$names=="otu"] = " "

########
###draw plots
##function that returns a ggplot for the given data frame (df)
##if legend is true, draw legend
getPlot <- function(df, legend) {
  df.m <- melt(df, id.vars="names") #columns to rows differentiated by measurement and value
  ylab = df$names
  plot <- ggplot(df.m, aes(x=variable, y=names)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low="white", high="darkgreen", limits=c(0, 4.5), na.value="lightgray", name="mean log\nrelative\nabundance") + 
    scale_x_discrete(name="", labels=c("rural\nT1", "rural\nT2", "urban \nT1 ", "urban\nT2")) +
    scale_y_discrete(name="", limits=rev(df$names), labels=rev(ylab)) +
    theme(axis.text.y=element_text(hjust=0, size=27, colour="black"), 
          axis.text.x=element_text(colour=c("blue", "blue", "red", "red"), size=24),
          # panel.background = element_rect(fill = 'grey', colour = 'white'),
          panel.grid.major = element_blank(),
          axis.line = element_blank(),
          axis.ticks.y = element_blank())
  if(!legend) { #first graph
    plot <- plot + annotate("text", label=taxaLevels, y=c(" ", "  ", "   "), x=2.5, size=15) +
      guides(fill=F)
  } else {
    plot <- plot + annotate("text", label="OTU", y=" ", x=2.5, size=15) +
      theme(legend.title=element_text(size=22), legend.text=element_text(size=22), legend.key.height=unit(30, "points"), legend.key.width=unit(20, "points"))
  }
  return(plot)
}

plot1 = getPlot(df1, F) #phylum, family, genus
plot2 = getPlot(df2, T) #OTU


###database figures
titleSize = 35
starSize = 65
xSize = 27
ySize = 24
yLabSize = 30
outlierSize = 3
boxSize = 1.5

##Silva
silva = read.table("silvaRuralVsUrban.txt",
                   header=T, colClasses = c(rep("numeric", 5), "logical"))
silva = silva[silva$pValue < 0.05,]
silva$ruralUrban = factor(ifelse(silva$higherInRural, "higher in\nrural", "higher in\nurban"))
s = ggplot(silva, aes(x=ruralUrban, y=percentIdentity)) + 
  stat_boxplot(geom ='errorbar', lwd=boxSize) + 
  geom_boxplot(outlier.size=outlierSize, lwd=boxSize) + #can use fatten to increase the median line
  labs(x="", y="Percent Identity") +
  ggtitle(expression(atop("SILVA database", atop("(1,756,783 sequences)")))) +
  theme(axis.text.x = element_text(size=xSize), axis.text.y = element_text(size=ySize),
        axis.title=element_text(size=yLabSize), plot.title = element_text(size=titleSize),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

##NCBI
ncbi = read.table("ncbiRuralVsUrban.txt",
                  header=T, colClasses = c("character", "character", rep("numeric", 5), "logical"))
ncbi = ncbi[ncbi$pValue < 0.05,]
ncbi$ruralUrban = factor(ifelse(ncbi$higherInRural, "higher in\nrural", "higher in\nurban"))
n = ggplot(ncbi, aes(x=ruralUrban, y=percentIdentity)) + 
  stat_boxplot(geom ='errorbar', lwd=boxSize) + 
  geom_boxplot(outlier.size=outlierSize, lwd=boxSize) + 
  labs(x="", y="Percent Identity") +
  ggtitle(expression(atop("NCBI Bacterial Genomes", atop("(2,767 sequences)")))) +
  theme(axis.text.x = element_text(size=xSize), axis.text.y = element_text(size=ySize),
        axis.title=element_text(size=yLabSize), plot.title = element_text(size=titleSize),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotation_custom(grob=textGrob(label="***", gp = gpar(fontsize = starSize)), ymin=103, ymax=99.99, xmin=1.5, xmax=1.5)
n = ggplotGrob(n) ##turn off clipping for stars
n$layout$clip[n$layout$name == "panel"] = "off"

##HMP
hmp = read.table("otusToMostWanted.txt",
                 header=T, sep="\t",
                 colClasses = c("character", "character", rep("numeric", 7), "character", "logical", "character"))
hmp = hmp[hmp$pValue < 0.05,]
hmp$ruralUrban = factor(ifelse(hmp$higherInRural, "higher in\nrural", "higher in\nurban"))
h1 = ggplot(hmp, aes(x=ruralUrban, y=percentIdentity)) + 
  stat_boxplot(geom ='errorbar', lwd=boxSize) + 
  geom_boxplot(outlier.size=outlierSize, lwd=boxSize) + 
  labs(x="", y="Percent Identity") +
  ggtitle("") +
  theme(axis.text.x = element_text(size=xSize), axis.text.y = element_text(size=ySize),
        axis.title=element_text(size=yLabSize),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
h2 = ggplot(hmp, aes(x=ruralUrban, y=hmp$stoolPrevelance)) + 
  stat_boxplot(geom ='errorbar', lwd=boxSize) + 
  geom_boxplot(outlier.size=outlierSize, lwd=boxSize) + 
  labs(x="", y="Gut Prevalence") +
  ggtitle("") +
  theme(axis.text.x = element_text(size=xSize), axis.text.y = element_text(size=ySize),
        axis.title=element_text(size=yLabSize), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotation_custom(grob=textGrob(label="****", gp = gpar(fontsize = starSize)), ymin=1.01, ymax=1.05, xmin=1.5, xmax=1.5)
h2 = ggplotGrob(h2)
h2$layout$clip[h2$layout$name == "panel"] = "off"

##combine plot
tiff("microbiomeHeatmapPlusDB.tif", res=350, height=7000, width=7000)
ggdraw() +
  draw_plot(plot1, x=0, y=.33, width=.345, height=.66) +
  draw_plot(plot2, x=.345, y=.33, width=.655, height=.66) +
  draw_plot(s, x=0, y=0, width=.25, height=.33) +
  draw_plot(n, x=.25, y=0, width=.25, height=.33) +
  draw_plot(h1, x=.5, y=0, width=.25, height=.29) +
  draw_plot(h2, x=.75, y=0, width=.25, height=.29) +
  draw_plot_label(label="Human Microbiome Project", x=.65, y=.31, size=titleSize, hjust=0, vjust=0, fontface="plain") +
  draw_plot_label(label=c("A", "B", "C", "D", "E"), size=45, hjust=0,
                  x=c(0, 0, .25, .5, .75), 
                  y=c(1, .307, .307, .307, .307))
dev.off()
