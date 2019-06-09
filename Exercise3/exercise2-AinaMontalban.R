#---------------------------------------------------------------------------------------------
##AINA MONTALBÁN
#---------------------------------------------------------------------------------------------


getwd()
#---------------------------------------------------------------------------------------------
###FOLDER DESTINATION DEFINITIONS
#---------------------------------------------------------------------------------------------
workingDir <-getwd()
dataDir <- file.path(workingDir, "dades")
resultsDir <- file.path(workingDir, "results")
setwd(resultsDir)


#---------------------------------------------------------------------------------------------
###INSTALLATION OF PACKAGES NEEDED
#---------------------------------------------------------------------------------------------
#if (!require(BiocManager)) install.packages("BiocManager")

installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
    BiocManager::install(pkg)
}else{
  require(pkg, character.only=T)
  }
}

#BiocManager::install("clariomsmousetranscriptcluster.db")

installifnot("pd.clariom.s.mouse")

installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("gridSVG")
#installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")


#---------------------------------------------------------------------------------------------
###LOAD DATA: TARGETS AND CEL FILES. 
#---------------------------------------------------------------------------------------------

#TARGETS
# Read the targets
targets <- read.csv(file = (file.path(dataDir,"targets.txt")),header = TRUE, sep = "")
targets

#CELFILES
# Read cel files
CELfiles <- list.celfiles(file.path(dataDir))
CELfiles
rawData <- read.celfiles(file.path(dataDir,CELfiles))

#DEFINE SOME VARIABLES FOR PLOTS
# Select sample names and colors
sampleNames <- as.character(targets$shortName)
sampleColor <- as.character(targets$Colors)
sampleNames
sampleColor
#---------------------------------------------------------------------------------------------
###QUALITY CONTROL OF ARRAYS: RAW DATA
#---------------------------------------------------------------------------------------------

#BOXPLOT
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)

#HIERARQUICAL CLUSTERING
#require(arrayQualityMetrics)
#arrayQualityMetrics(rawData, outdir = file.path("./results", 272
                                                #"QCDir.Raw"), force=TRUE)


clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)

#PRINCIPAL COMPONENT ANALYSIS
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

#SAVE TO A FILE
pdf(file.path(resultsDir, "QCPlots_Raw.pdf"))
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples of RawData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()


#---------------------------------------------------------------------------------------------
###DATA NORMALIZATION
#---------------------------------------------------------------------------------------------
eset_rma<-rma(rawData)

write.exprs(eset_rma, file.path(resultsDir, "NormData.txt"))


  #---------------------------------------------------------------------------------------------
###QUALITY CONTROL OF ARRAYS: NORMALIZED DATA
#---------------------------------------------------------------------------------------------

#BOXPLOT
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)

#HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(eset))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)

#PRINCIPAL COMPONENT ANALYSIS
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(eset), labels=sampleNames, dataDesc="NormData", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

#SAVE TO A FILE
pdf(file.path(resultsDir, "QCPlots_Norm.pdf"))
boxplot(eset_rma, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(eset_rma), labels=sampleNames, dataDesc="selected samples", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()

#ARRAY QUALITY METRICS
arrayQualityMetrics(eset,  reporttitle="QualityControl", force=TRUE)






#---------------------------------------------------------------------------------------------
###Detecting Most Variable Genes
#---------------------------------------------------------------------------------------------

sds <- apply (exprs(eset_rma), 1, sd) 
sdsO<- sort(sds) 
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",sub="Vertical lines represent 90% and 95% percentiles",
       xlab="Gene index (from least to most variable)",
       ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))

#---------------------------------------------------------------------------------------------
###FILTER OUT THE DATA
#---------------------------------------------------------------------------------------------

## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 29129 features, 20 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: GSM3559494 GSM3559495 ... GSM3559513 (20 total)
##   varLabels: sampleName genotype replicate
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation: GPL23038
gse <- getGEO("GSE124927")
eset <- gse[[1]]
annotation(eset) <- "clariomsmousetranscriptcluster.db"
eset_filtered <- nsFilter(eset, var.func=IQR,
         var.cutoff=0.75, var.filter=TRUE,
         filterByQuantile=TRUE)
#NUMBER OF GENES OUT
print(eset_filtered$filter.log$numLowVar)
print(eset_filtered$filter.log$numDupsRemoved)

#29129 - eset_filtered$filter.log$numLowVar
#NUMBER OF GENES IN
print(eset_filtered$eset)

#---------------------------------------------------------------------------------------------
###DIFERENTIAL EXPRESSED GENES SELECTION. LINEAR MODELS. COMPARITIONS
#---------------------------------------------------------------------------------------------

#CONTRAST MATRIX.lINEAR MODEL
treat <- targets$shortName
treat
lev <- factor(treat, levels = unique(treat))
lev
design <-model.matrix(~0+lev)
design
colnames(design) <- levels(lev)
rownames(design) <- sampleNames
print(design)

#COMPARISON
cont.matrix1 <- makeContrasts( 
  CD.vs.HF = CD- HF, 
  HF.vs.HF_RES = HF-HF_RES, 
  CD.vs.HF_RES = CD-HF_RES,
  levels = design)
cont.matrix1


cont.matrix2 <- makeContrasts( 
  CD.vs.HF = HF- CD, 
  HF.vs.HF_RES = HF_RES-HF, 
  CD.vs.HF_RES = HF_RES-CD,
  levels = design)
cont.matrix2
#comparison1 <- "Effect of Diet"

#MODEL FIT
fit1 <- lmFit(eset_filtered$eset, design)
fit.main1 <- contrasts.fit(fit1, cont.matrix1)
fit.main1 <- eBayes(fit.main1)

#---------------------------------------------------------------------------------------------
###DIFERENTIAL EXPRESSED GENES LISTS.TOPTABLES
#---------------------------------------------------------------------------------------------
comparison1 <- "CDvsHF"
#FILTER BY FALSE DISCOVERY RATE AND FOLD CHANGE
topTab_CDvsHF <- topTable(fit.main1, number=nrow(fit.main1), coef = "CD.vs.HF", adjust="fdr")
head(topTab_CDvsHF)

comparison2 <- "HFvsHF_RES"
#minimum absolute log2-fold-change required. 
topTab_HFvsHF_RES <- topTable(fit.main1, number=nrow(fit.main1), coef = "HF.vs.HF_RES", adjust="fdr")
head(topTab_HFvsHF_RES)

comparison3 <- "CDvsHF_RES"
topTab_CDvsHF_RES <- topTable(fit.main1, number=nrow(fit.main1), coef = "CD.vs.HF_RES", adjust="fdr")
head(topTab_CDvsHF_RES)


#EXPORTED TO CSV AND HTML FILE
write.csv2(topTab_CDvsHF, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison1, ".csv", sep = "")))

write.csv2(topTab_HFvsHF_RES, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison2, ".csv", sep = "")))

write.csv2(topTab_CDvsHF_RES, file= file.path(resultsDir,paste("Selected.Genes.in.comparison.",
                                                    comparison3, ".csv", sep = "")))


print(xtable(topTab_CDvsHF,align="lllllll"),type="html",html.table.attributes="",
      file=paste("Selected.Genes.in.comparison.",comparison1,".html", sep=""))

print(xtable(topTab_HFvsHF_RES,align="lllllll"),type="html",html.table.attributes="",
      file=paste("Selected.Genes.in.comparison.",comparison2,".html", sep=""))

print(xtable(topTab_CDvsHF_RES,align="lllllll"),type="html",html.table.attributes="",
      file=paste("Selected.Genes.in.comparison.",comparison3,".html", sep=""))

#---------------------------------------------------------------------------------------------
###VOLCANO PLOTS
#---------------------------------------------------------------------------------------------
volcanoplot(fit.main1, highlight=10, names=fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep="\n"))
abline(v = c(-3, 3))


pdf(file.path(resultsDir,"Volcanos.pdf"))
volcanoplot(fit.main1, highlight = 10, names = fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep = "\n"))
abline(v = c(-3, 3))
dev.off()


#---------------------------------------------------------------------------------------------
###HEATMAP PLOTS
#---------------------------------------------------------------------------------------------
topTab.heatmap1 <- topTab.end1
rownames(topTab.heatmap1) <- topTab.heatmap1[,1]
topTab.heatmap1_subset <- topTab.heatmap1[,c(5:24)]

#PREPARE THE DATA
my_frame <- data.frame(exprs(eset_rma))
head(my_frame)
HMdata <- merge(my_frame, topTab.heatmap1_subset, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1, 2:20)]
head(HMdata)
HMdata2 <- data.matrix(HMdata, rownames.force=TRUE)
head(HMdata2)
write.csv2(HMdata2, file = file.path(resultsDir,"Data2HM.csv"))
traceback()
#HEATMAP PLOT
my_palette <- colorRampPalette(c("blue", "red","green"))(n = 299)


#EXPORT TO PDF FILE
pdf(file.path(resultsDir,"HeatMap_CDvsHF_RES.pdf"))
heatmap.2(HMdata2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap_CDvsHF_RES FC>=3",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=c("red", "blue", "red", "blue", 
                          "green", "green", "blue", "blue",
                          "red", "red", "red", "red","blue", 
                          "blue", "red", "green", "green",
                          "red","red","red"), 
          tracecol=NULL,
          srtCol=30)
dev.off()

topTab.heatmap2 <- topTab.end2
rownames(topTab.heatmap2) <- topTab.heatmap2[,1]
topTab.heatmap2_subset <- topTab.heatmap2[,c(5:24)]

#PREPARE THE DATA
my_frame <- data.frame(exprs(eset_rma))
head(my_frame)
HMdata_2 <- merge(my_frame, topTab.heatmap2_subset, by.x = 0, by.y = 0)
rownames(HMdata_2) <- HMdata_2$Row.names
HMdata_2 <- HMdata_2[, -c(1, 2:20)]
head(HMdata_2)
HMdata2_2 <- data.matrix(HMdata_2, rownames.force=TRUE)
head(HMdata2_2)
write.csv2(HMdata2_2, file = file.path(resultsDir,"Data2HM_2.csv"))
traceback()
#HEATMAP PLOT
my_palette <- colorRampPalette(c("blue", "red","green"))(n = 299)

head(exprs(rawData))
#EXPORT TO PDF FILE
pdf(file.path(resultsDir,"HeatMap_CDvsHF_RES.pdf"))
heatmap.2(HMdata2_2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap_CDvsHF_RES FC>=3",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=c("red", "blue", "red", "blue", 
                          "green", "green", "blue", "blue",
                          "red", "red", "red", "red","blue", 
                          "blue", "red", "green", "green",
                          "red","red","red"), 
          tracecol=NULL,
          srtCol=30)
dev.off()








traceback()
#---------------------------------------------------------------------------------------------
###DATA ANNOTATION

#-----------------------------------------------------
----------------------------------------
head(ls("package:clariomsmousetranscriptcluster.db"))
head(ls("package:org.Mm.eg.db"))
# > head(ls("package:clariomsmousetranscriptcluster.db"))
#[1] "clariomsmousetranscriptcluster"            "clariomsmousetranscriptclusterACCNUM"     
#[3] "clariomsmousetranscriptclusterALIAS2PROBE" "clariomsmousetranscriptclusterCHR"        
#[5] "clariomsmousetranscriptclusterCHRLENGTHS"  "clariomsmousetranscriptclusterCHRLOC" 

### ALL THE GENES
all_anota<-data.frame(exprs(eset_rma))
Annot <- data.frame(SYMBOL=sapply(contents(clariomsmousetranscriptclusterSYMBOL), paste, collapse=", "),
                    DESC=sapply(contents(clariomsmousetranscriptclusterGENENAME), paste, collapse=", "),
                    UNIPROT=sapply(contents(clariomsmousetranscriptclusterUNIPROT), paste, collapse=", "))


Annot<-Annot[!Annot$SYMBOL=="NA",]
Annot<-Annot[!Annot$DESC=="NA",]
Annot<-Annot[!Annot$UNIPROT=="NA",]

head(Annot)

anotaGenes <- merge(Annot,all_anota, by.x=0,by.y=0)
head(anotaGenes)
write.table(anotaGenes, file ="data.ann.txt",sep="\t")

rownames(anotaGenes) <- anotaGenes[,1]
anotaGenes <- anotaGenes[,-1]


topTab.end1 <- anotaGenes.end1[,c(1:4,5:24)]
rownames(topTab.end1) <- topTab.end1[,1]

## Comparison 1
anotaGenes.end1 <- merge(anotaGenes, topTab_CDvsHF, by.x=0,by.y=0)
topTab.end1 <- merge(anotaGenes, topTab_CDvsHF, by.x=0,by.y=0)
#reordenamos las columnas
#topTab.end1 <- anotaGenes.end1[,c(1:3,12:17,4:11)]
topTab.end1 <- topTab.end1[order(-topTab.end1$B),]

rownames(topTab.end1) <- topTab.end1[,1]
topTab.end1 <- topTab.end1[, -1]
write.csv(topTab.end1, file = file.path(resultsDir,"TopTable_comp1.end.csv"))



## Comparison 2
anotaGenes.end2 <- merge(anotaGenes, topTab_HFvsHF_RES, by.x=0,by.y=0)
#reordenamos las columnas
topTab.end2 <- anotaGenes.end2
topTab.end2 <- topTab.end2[order(-topTab.end2$B),]

rownames(topTab.end2) <- topTab.end2$Row.names
#topTab.end2 <- topTab.end2[, -1]
write.csv(topTab.end2, file = file.path(resultsDir,"TopTable_comp2.end.csv"))

## Comparison 3
anotaGenes.end3 <- merge(anotaGenes, topTab_HFvsHF_RES, by.x=0,by.y=0)
#reordenamos las columnas
topTab.end3 <- anotaGenes.end3
topTab.end3 <- topTab.end3[order(-topTab.end3$B),]

rownames(topTab.end3) <- topTab.end3[,1]
topTab.end3 <- topTab.end3[, -1]
write.csv(topTab.end3, file = file.path(resultsDir,"TopTable_comp3.end.csv"))


#---------------------------------------------------------------------------------------------
#END OF SCRIPT
#---------------------------------------------------------------------------------------------
#######

annotatedTopTable <- function(topTab, anotPackage){
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID","GENENAME", "UNIPROT"))  
  annotatedTopTab <- merge(x=geneAnots, y=topTab, by.x = "PROBEID",by.y="PROBEID")
  return(annotatedTopTab)
}

topAnnotated_CDvsHF <- annotatedTopTable(topTab_CDvsHF, "clariomsmousetranscriptcluster.db")
topAnnotated_HFvsHF_RES <- annotatedTopTable(topTab_HFvsHF_RES, "clariomsmousetranscriptcluster.db")
topAnnotated_CDvsHF_RES <- annotatedTopTable(topTab_CDvsHF_RES, "clariomsmousetranscriptcluster.db")

########




geneSym <- select(clariomsmousetranscriptcluster.db, rownames(fit.main1), c("SYMBOL"))  
SYMBOLS <- geneSym$SYMBOL

cont.matrix1[1]
pdf(file.path(resultsDir,"Volcanos_Groups.pdf"))
volcanoplot(fit.main1, coef = 1, highlight = 10, names = SYMBOLS, 
            main = paste("Differentially expressed genes", 
            colnames(cont.matrix1)[1], sep = "\n"))
abline(v = c(-1, 1))
volcanoplot(fit.main1, coef = 2, highlight = 10, names = SYMBOLS, 
            main = paste("Differentially expressed genes", 
                         colnames(cont.matrix1)[2], sep = "\n"))
abline(v = c(-1, 1))
volcanoplot(fit.main1, coef =3, highlight = 10, names = SYMBOLS, 
            main = paste("Differentially expressed genes", 
                         colnames(cont.matrix1)[3], sep = "\n"))
abline(v = c(-1, 1))
dev.off()



### Multiple Comparisons
fit.main1$p.value
res <- decideTests(fit.main1, method = "separate", adjust.method = "fdr", p.value = 0.7)
sum.res.rows <- apply(abs(res), 1, sum)
sum.res.rows
res.selected <- res[sum.res.rows!=0,]
res.selected
print(summary(res))


