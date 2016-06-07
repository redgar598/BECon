##################### Get cell composition probes
setwd("~/Documents/Blood_Brain/")
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(mclust)
library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kmanifest)
library(quadprog)
load("BLBR_Beta_Combatted_Mval_bmiqonly.RData") #441198     70
          RGset = as.data.frame(combat_BLBR_Beta)[,1:3]



#meanPlot = TRUE
cellType = "Blood"
verbose = TRUE

  require(quadprog)
  library(genefilter)
  library(matrixStats)
  platform <- "450k"
  referencePkg <- sprintf("FlowSorted.%s.%s", cellType, platform)
  if (!require(referencePkg, character.only = TRUE)) 
    stop(sprintf("Could not find reference data package for cellType '%s' and platform '%s' (inferred package name is '%s')", 
                 cellType, platform, referencePkg))
  data(list = referencePkg)
  referenceRGset <- get(referencePkg)
  if (verbose) 
    cat("[estimateCellCounts] Combining Data with Flow Sorted Data.\n")
refpd <- data.frame(sampleNames = sampleNames(referenceRGset), 
                    studyIndex = rep("reference", times = ncol(referenceRGset)), 
                    stringsAsFactors = FALSE)
refpd <- cbind(refpd,referenceRGset$CellType)
referenceMset<-as.data.frame(getBeta(referenceRGset))
compData <- pickCompProbes(referenceMset,refpd)
               

pickCompProbes <- function(Mset, pdref, numProbes = 50) {
  splitit <- function(x) {
    split(seq(along=x), x)
  }
  
  #p <- getBeta(Mset)
  p <- as.matrix(Mset)
  #pd <- as.data.frame(pData(Mset))
  pd <- as.data.frame(pdref)
  
  ## only keep 6 components from kere
  keep <- which(pd[,3] %in% c("Mono", "Bcell", 
                              "Gran", "CD4T", "CD8T", "NK"))
  pd <- pd[keep,]
  p <- p[,keep]
  
  ## make cell type a factor 
  pd$CellType <- factor(pd[,3], 
                        levels = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"))
  
  # get fstats
  ffComp <- rowFtests(p, pd$CellType)
  prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
  r <- rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
  names(compTable)[c(1, 9:11)] <- c("Fstat", "low", "high", "range")
  
  # t-test by cell type
  tIndexes <- splitit(pd$CellType)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  
  ## take N up and N down
  probeList <- lapply(tstatList, function(x) {
    y <- x[x[,"p.value"] < 1e-8,]
    yUp <- y[order(y[,"dm"], decreasing=TRUE),]
    yDown <- y[order(y[,"dm"], decreasing=FALSE),]
    c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
  })
  
  trainingProbes <- unlist(probeList)
  p <- p[trainingProbes,]
  
  pMeans <- colMeans(p)
  names(pMeans) <- pd$CellType
  
  mod <- model.matrix(~pd$CellType-1)
  colnames(mod) <- levels(pd$CellType)
  #form <- as.formula(sprintf("y ~ %s - 1", colnames(mod), collapse="+"))
  form <- as.formula("y ~ CD8T + CD4T + NK + Bcell + Mono + Gran -1")
  
  tmp <- validationWBC(p,data.frame(mod),form)
  coefEsts <- tmp$coefEsts
  
  out <- list(coefEsts = coefEsts, compTable = compTable,
              sampleMeans = pMeans)
  return(out)
}
Blood_training_probes<-rownames(coefs)

save(Blood_training_probes, file="Blood_training_probes.RData")



##### Genes associated
CETS<-read.csv("CETS.sites.csv")
CETS_probes<-as.character(CETS$x)


load("Gene_CpG_Relations_updatejune2015.RData")
Gene_CpG_Relations_update<-Gene_CpG_Relations_update[!duplicated(Gene_CpG_Relations_update), ]

Genes_correlated_CpGs_CETS<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%CETS_probes),]
Genes_correlated_CpGs_CETS<-Genes_correlated_CpGs_CETS[,c(1:4,6)]
Genes_correlated_CpGs_CETS<-Genes_correlated_CpGs_CETS[!duplicated(Genes_correlated_CpGs_CETS),]
Genes_correlated_CpGs_CETS<-Genes_correlated_CpGs_CETS[!duplicated(Genes_correlated_CpGs_CETS[,c(1,4)]),]#remove duplicate CpG to gene associations
length(unique(Genes_correlated_CpGs_CETS$Probe_ID))#10000 CpGs
length(unique(Genes_correlated_CpGs_CETS$gene))#4972 Genes

Genes_correlated_CpGs_blood<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%Blood_training_probes),]
Genes_correlated_CpGs_blood<-Genes_correlated_CpGs_blood[,c(1:4,6)]
Genes_correlated_CpGs_blood<-Genes_correlated_CpGs_blood[!duplicated(Genes_correlated_CpGs_blood),]
Genes_correlated_CpGs_blood<-Genes_correlated_CpGs_blood[!duplicated(Genes_correlated_CpGs_blood[,c(1,4)]),]#remove duplicate CpG to gene associations
length(unique(Genes_correlated_CpGs_blood$Probe_ID))#600 CpGs
length(unique(Genes_correlated_CpGs_blood$gene))#426 Genes


write.csv(Genes_correlated_CpGs_CETS, file="Methylation_Brain_CETS_deconvolution_sites.csv")
write.csv(Genes_correlated_CpGs_blood, file="Methylation_blood_Houseman_deconvolution_sites.csv")