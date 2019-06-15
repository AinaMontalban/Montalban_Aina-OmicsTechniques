## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------
workingDir <-getwd()
dataDir <- file.path(workingDir, "dades")
resultsDir <- file.path(workingDir, "results")
setwd(resultsDir)

## ----message=FALSE-------------------------------------------------------
#if (!require(BiocManager)) install.packages("BiocManager")

installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
    BiocManager::install(pkg)
}else{
  require(pkg, character.only=T)
  }
}

installifnot("pd.clariom.s.mouse")
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("gridSVG")
installifnot("ReactomePA")
#installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")



## ----Data Capture TARGETS------------------------------------------------
targets <- read.csv(file = (file.path(dataDir,"targets.txt")),header = TRUE, sep = "")
targets


## ----Data Capture CEL, message=FALSE, results='hide'---------------------
CELfiles <- list.celfiles(file.path(dataDir))
CELfiles
rawData <- read.celfiles(file.path(dataDir,CELfiles))


## ------------------------------------------------------------------------
# Select sample names and colors
sampleNames <- as.character(targets$shortName)
sampleColor <- as.character(targets$Colors)


## ------------------------------------------------------------------------
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
       cex.axis=0.6, col=sampleColor, names=sampleNames)


## ------------------------------------------------------------------------
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
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
        formapunts=c(rep(16,4)), myCex=0.6)




## ------------------------------------------------------------------------
eset_rma<-rma(rawData)
write.exprs(eset_rma, file.path(resultsDir, "NormData.txt"))


## ------------------------------------------------------------------------
boxplot(eset_rma, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)


## ------------------------------------------------------------------------
clust.euclid.average <- hclust(dist(t(exprs(eset_rma))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
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

plotPCA(exprs(eset_rma), labels=sampleNames, dataDesc="NormData", colors=sampleColor,
        formapunts=c(rep(16,4)), myCex=0.6)


## ------------------------------------------------------------------------
sds <- apply (exprs(eset_rma), 1, sd) 
sdsO<- sort(sds) 
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",sub="Vertical lines represent 90% and 95% percentiles",
       xlab="Gene index (from least to most variable)",
       ylab="Standard deviation")
#abline(v=length(sds)*c(0.9,0.95))


## ------------------------------------------------------------------------
annotation(eset_rma) <- "clariomsmousetranscriptcluster.db"
eset_filtered <- nsFilter(eset_rma, var.func=IQR,
         var.cutoff=0.75, var.filter=TRUE,
         filterByQuantile=TRUE)
print(eset_filtered$filter.log$numLowVar)
print(eset_filtered$eset)


## ------------------------------------------------------------------------
#CONTRAST MATRIX.lINEAR MODEL
treat <- targets$shortName
lev <- factor(treat, levels = unique(treat))
design <-model.matrix(~ 0+targets$shortName)
rownames(design)<- targets$sampleName 
colnames(design)<-c("CD", "HF", "HF_RES")
print(design)


## ------------------------------------------------------------------------
cont.matrix1 <- makeContrasts( 
  CD.vs.HF = CD- HF, 
  HF.vs.HF_RES = HF-HF_RES, 
  CD.vs.HF_RES = CD-HF_RES,
  levels = design)
cont.matrix1


## ------------------------------------------------------------------------
#MODEL FIT
fit1 <- lmFit(eset_filtered$eset, design)
fit.main1 <- contrasts.fit(fit1, cont.matrix1)
fit.main1 <- eBayes(fit.main1)


## ------------------------------------------------------------------------
comparison1 <- "CDvsHF"
#FILTER BY FALSE DISCOVERY RATE AND FOLD CHANGE
topTab_CDvsHF <- topTable(fit.main1, number=nrow(fit.main1), coef = "CD.vs.HF", adjust="fdr", lfc = 1)


## ------------------------------------------------------------------------
head(topTab_CDvsHF)


## ------------------------------------------------------------------------
comparison2 <- "HFvsHF_RES"
#minimum absolute log2-fold-change required. 
topTab_HFvsHF_RES <- topTable(fit.main1, number=nrow(fit.main1), coef = "HF.vs.HF_RES", adjust="fdr", lfc = 1)


## ------------------------------------------------------------------------
head(topTab_HFvsHF_RES)


## ------------------------------------------------------------------------
comparison3 <- "CDvsHF_RES"
topTab_CDvsHF_RES <- topTable(fit.main1, number=nrow(fit.main1), coef = "CD.vs.HF_RES", adjust="fdr", lfc = 1)


## ------------------------------------------------------------------------
head(topTab_CDvsHF_RES)


## ------------------------------------------------------------------------
all_anota<-data.frame(exprs(eset_rma))
Annot <- data.frame(SYMBOL=sapply(contents(clariomsmousetranscriptclusterSYMBOL), paste, collapse=", "),
                    DESC=sapply(contents(clariomsmousetranscriptclusterGENENAME), paste, collapse=", "),
                    UNIPROT=sapply(contents(clariomsmousetranscriptclusterUNIPROT), paste, collapse=", "),
                    ENTREZID=sapply(contents(clariomsmousetranscriptclusterENTREZID), paste, collapse=", "))
Annot<-Annot[!Annot$SYMBOL=="NA",]
Annot<-Annot[!Annot$DESC=="NA",]
Annot<-Annot[!Annot$UNIPROT=="NA",]
Annot<-Annot[!Annot$ENTREZID=="NA",]
anotaGenes <- merge(Annot,all_anota, by.x=0,by.y=0)
write.table(anotaGenes, file ="data.ann.txt",sep="\t")
rownames(anotaGenes) <- anotaGenes[,1]
anotaGenes <- anotaGenes[,-1]


## ------------------------------------------------------------------------
## Comparison 1
anotaGenes.end1 <- merge(anotaGenes, topTab_CDvsHF, by.x=0,by.y=0)
topTab.end1 <- merge(anotaGenes, topTab_CDvsHF, by.x=0,by.y=0)
#reordenamos las columnas
#topTab.end1 <- anotaGenes.end1[,c(1:3,12:17,4:11)]
topTab.end1 <- topTab.end1[order(-topTab.end1$B),]
rownames(topTab.end1) <- topTab.end1[,1]
topTab.end1 <- topTab.end1[, -1]
write.csv(topTab.end1, file = file.path(resultsDir,"TopTable_comp1.end.csv"))


## ------------------------------------------------------------------------
head(topTab.end1[,1:3])


## ------------------------------------------------------------------------
## Comparison 2
anotaGenes.end2 <- merge(anotaGenes, topTab_HFvsHF_RES, by.x=0,by.y=0)
#reordenamos las columnas
topTab.end2 <- anotaGenes.end2
topTab.end2 <- topTab.end2[order(-topTab.end2$B),]

rownames(topTab.end2) <- topTab.end2[,1]
topTab.end2 <- topTab.end2[, -1]
write.csv(topTab.end2, file = file.path(resultsDir,"TopTable_comp2.end.csv"))


## ------------------------------------------------------------------------
head(topTab.end2[,1:3])


## ------------------------------------------------------------------------
## Comparison 3
anotaGenes.end3 <- merge(anotaGenes, topTab_HFvsHF_RES, by.x=0,by.y=0)
#reordenamos las columnas
topTab.end3 <- anotaGenes.end3
topTab.end3 <- topTab.end3[order(-topTab.end3$B),]

rownames(topTab.end3) <- topTab.end3[,1]
topTab.end3 <- topTab.end3[, -1]
write.csv(topTab.end3, file = file.path(resultsDir,"TopTable_comp3.end.csv"))


## ------------------------------------------------------------------------
head(topTab.end3[,1:3])


## ------------------------------------------------------------------------
geneSym <- select(clariomsmousetranscriptcluster.db, rownames(fit.main1), c("SYMBOL"))  
SYMBOLS <- geneSym$SYMBOL


volcanoplot(fit.main1, coef = 1, highlight = 10, names = SYMBOLS, 
            main = paste("Differentially expressed genes", 
            colnames(cont.matrix1)[1], sep = "\n"))
abline(v = c(-1, 1))
volcanoplot(fit.main1, coef = 2, highlight = 10, names = SYMBOLS, 
            main = paste("Differentially expressed genes", 
                         colnames(cont.matrix1)[2], sep = "\n"))
abline(v = c(-1, 1))
volcanoplot(fit.main1, coef = 3, highlight = 10, names = SYMBOLS, 
            main = paste("Differentially expressed genes", 
                         colnames(cont.matrix1)[3], sep = "\n"))
abline(v = c(-1, 1))


## ------------------------------------------------------------------------
res <- decideTests(fit.main1, method = "separate", adjust.method = "fdr", p.value = 1, lfc = 1)
sum.res.rows <- apply(abs(res), 1, sum)
res.selected <- res[sum.res.rows!=0,]
print(summary(res))
vennDiagram(res.selected[,1:3], cex=0.9)

