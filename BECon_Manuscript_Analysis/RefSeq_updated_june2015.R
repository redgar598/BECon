##############
## CpG Gene Features (with ncRNA)
##############

setwd("/big_data/redgar/BECon")
library(ggplot2)
library(IlluminaHumanMethylation450k.db)
library(reshape)
library(RColorBrewer)

    # Data
    ucsc<-read.table("Gene_annotation/refseq") # 46722 gene bodies
    ucsc<-ucsc[,c(3,4,5,6,13)]
    colnames(ucsc)<-c("chr","strand","txStart","txEnd", "name2")
    refseq<-unique(ucsc)#33431 genes 

# Merge Overlaping refseq genes, Come back to Uusc file for acessions
refseq<-split(ucsc, ucsc$name2)
refseq_no_duplicates<-lapply(refseq, function(x) x[!duplicated(x),])
refseq_isoforms<-lapply(refseq_no_duplicates, function(x) {x$isoform<-seq(1, nrow(x), 1)
                                                           x})

refseq<-do.call(rbind, refseq_isoforms)#24047 genes; isoforms 33431
save(refseq, file="Gene_annotation/refseq_with_isoforms.RData")

# Strand Issues, txstart is just lowest number coordinate, in minus strand genes txstart is actually gene end
# Flip the minus strand coordiantes
refseq_plus<-subset(refseq, strand=="+")
refseq_minus<-subset(refseq, strand=="-")
refseq_minus$fixStart<-refseq_minus$txEnd
refseq_minus$txEnd<-refseq_minus$txStart
refseq_minus$txStart<-refseq_minus$fixStart
refseq_minus$fixStart<-NULL

# Chromosome Issues
# CpGI 
chr <- IlluminaHumanMethylation450kCHR37
chr <- as.data.frame(chr)

coor <- IlluminaHumanMethylation450kCPGCOORDINATE
mapped_probes <- mappedkeys(coor)
coor<- as.data.frame(coor[mapped_probes])

CpG<-merge(coor, chr, by="Probe_ID")
CpG<-CpG[which(CpG$Chromosome_37!=""),] ## Remove CpGs with no chromosome leaving 485512 CpGs

CpG$Chromosome_37<-as.factor(CpG$Chromosome_37)
levels(CpG$Chromosome_37)<-c(levels(CpG$Chromosome_37)[1:22], 23,24)
CpG$Chromosome_37<-as.factor(as.character(CpG$Chromosome_37))
CpG<-CpG[order(CpG$Chromosome_37),]
CpG<-split(CpG, CpG$Chromosome_37)
# RefSeq
#in this update I have removed all haplotype chromosomes since they cause craziness
levels(refseq_plus$chr)<-c(1,10,11,12,13,14,15,16,17,0,0,18,19,0,0,0,2,20,21,22,3,4,0,0,0,5,6,0,0,0,0,0,0,0,7,0,8,9,0,0,0,0,0,0,0,0,0,23,24)
refseq_plus<-subset(refseq_plus, chr!=0)
refseq_plus$chr<-as.factor(as.character(refseq_plus$chr))
refseq_plus<-refseq_plus[order(refseq_plus$chr),]
refseq_plus<-split(refseq_plus, refseq_plus$chr)
levels(refseq_minus$chr)<-c(1,10,11,12,13,14,15,16,17,0,0,18,19,0,0,0,2,20,21,22,3,4,0,0,0,5,6,0,0,0,0,0,0,0,7,0,8,9,0,0,0,0,0,0,0,0,0,23,24)
refseq_minus<-subset(refseq_minus, chr!=0)
refseq_minus$chr<-as.factor(as.character(refseq_minus$chr))
refseq_minus<-refseq_minus[order(refseq_minus$chr),]
refseq_minus<-split(refseq_minus, refseq_minus$chr)

##############
## Catagorize CpGs
##############

# CpG is promoter from 300bp before the tss to 1500 bp after, 
#intragenic from 300bp after tss to 300 bp before the txend, 
#3' if a CpG is 300bp before the txend or 300bp after,
#intergenic if a CpG is in none of these catagories on any gene


# Promoter CpGs
promoter_plus<-lapply(1:length(CpG),function(chr) do.call(rbind,lapply(1:nrow(refseq_plus[[chr]]), function(y) {
  probes_on_gene<-CpG[[chr]][which(CpG[[chr]]$Coordinate_37<=refseq_plus[[chr]]$txStart[y]+300 & 
                                     CpG[[chr]]$Coordinate_37>=refseq_plus[[chr]]$txStart[y]-1500),]
  if(nrow(probes_on_gene)>0){probes_on_gene$gene<-refseq_plus[[chr]]$name2[y]
                             probes_on_gene$isoform<-refseq_plus[[chr]]$isoform[y]}else{}
  if(nrow(probes_on_gene)>0){probes_on_gene$region<-"promoter"}else{}
  probes_on_gene})))

promoter_minus<-lapply(1:length(CpG),function(chr) do.call(rbind,lapply(1:nrow(refseq_minus[[chr]]), function(y) {
  probes_on_gene<-CpG[[chr]][which(CpG[[chr]]$Coordinate_37>=refseq_minus[[chr]]$txStart[y]-300 & 
                                     CpG[[chr]]$Coordinate_37<=refseq_minus[[chr]]$txStart[y]+1500),]
  if(nrow(probes_on_gene)>0){probes_on_gene$gene<-refseq_minus[[chr]]$name2[y]
                             probes_on_gene$isoform<-refseq_minus[[chr]]$isoform[y]}else{}
  if(nrow(probes_on_gene)>0){probes_on_gene$region<-"promoter"}else{}
  probes_on_gene})))


promoter_plus<-do.call(rbind, promoter_plus)
promoter_minus<-do.call(rbind, promoter_minus)
promoter<-rbind(promoter_plus,promoter_minus) #47.5% 230580


# Intragenic
intragenic_plus<-lapply(1:length(CpG),function(chr) do.call(rbind,lapply(1:nrow(refseq_plus[[chr]]), function(y) {
  probes_on_gene<-CpG[[chr]][which(CpG[[chr]]$Coordinate_37>=refseq_plus[[chr]]$txStart[y]+300 & 
                                     CpG[[chr]]$Coordinate_37<=refseq_plus[[chr]]$txEnd[y]-300),]
  if(nrow(probes_on_gene)>0){probes_on_gene$gene<-refseq_plus[[chr]]$name2[y]
                             probes_on_gene$isoform<-refseq_plus[[chr]]$isoform[y]}else{}
  if(nrow(probes_on_gene)>0){probes_on_gene$region<-"intragenic"}else{}
  probes_on_gene})))

intragenic_minus<-lapply(1:length(CpG),function(chr) do.call(rbind,lapply(1:nrow(refseq_minus[[chr]]), function(y) {
  probes_on_gene<-CpG[[chr]][which(CpG[[chr]]$Coordinate_37<=refseq_minus[[chr]]$txStart[y]-300 & 
                                     CpG[[chr]]$Coordinate_37>=refseq_minus[[chr]]$txEnd[y]+300),]
  if(nrow(probes_on_gene)>0){probes_on_gene$gene<-refseq_minus[[chr]]$name2[y]
                             probes_on_gene$isoform<-refseq_minus[[chr]]$isoform[y]}else{}
  if(nrow(probes_on_gene)>0){probes_on_gene$region<-"intragenic"}else{}
  probes_on_gene})))

intragenic_plus<-do.call(rbind, intragenic_plus)
intragenic_minus<-do.call(rbind, intragenic_minus)
intragenic<-rbind(intragenic_plus,intragenic_minus) #68.0% 330214

# 3'
three_plus<-lapply(1:length(CpG),function(chr) do.call(rbind,lapply(1:nrow(refseq_plus[[chr]]), function(y) {
  probes_on_gene<-CpG[[chr]][which(CpG[[chr]]$Coordinate_37>=refseq_plus[[chr]]$txEnd[y]-300 & 
                                     CpG[[chr]]$Coordinate_37<=refseq_plus[[chr]]$txEnd[y]+300),]
  if(nrow(probes_on_gene)>0){probes_on_gene$gene<-refseq_plus[[chr]]$name2[y]
                             probes_on_gene$isoform<-refseq_plus[[chr]]$isoform[y]}else{}
  if(nrow(probes_on_gene)>0){probes_on_gene$region<-"three_plus"}else{}
  probes_on_gene})))

three_minus<-lapply(1:length(CpG),function(chr) do.call(rbind,lapply(1:nrow(refseq_minus[[chr]]), function(y) {
  probes_on_gene<-CpG[[chr]][which(CpG[[chr]]$Coordinate_37<=refseq_minus[[chr]]$txEnd[y]+300 & 
                                     CpG[[chr]]$Coordinate_37>=refseq_minus[[chr]]$txEnd[y]-300),]
  if(nrow(probes_on_gene)>0){probes_on_gene$gene<-refseq_minus[[chr]]$name2[y]
                             probes_on_gene$isoform<-refseq_minus[[chr]]$isoform[y]}else{}
  if(nrow(probes_on_gene)>0){probes_on_gene$region<-"three_plus"}else{}
  probes_on_gene})))

three_plus<-do.call(rbind, three_plus)
three_minus<-do.call(rbind, three_minus)
three_prime<-rbind(three_plus,three_minus) #3.65% 16762

CpG_all<-do.call(rbind, CpG)
intergenic<-CpG_all[which(!(CpG_all$Probe_ID%in%c(three_plus$Probe_ID, intragenic$Probe_ID, promoter$Probe_ID))),] 
intergenic$gene<-"None"
intergenic$isoform<-"None"
intergenic$region<-"intergenic" #24.3% 111353

Gene_CpG_Relations_update<-rbind(three_prime, intragenic,promoter, intergenic) # 688909 CpG gene relations
save(Gene_CpG_Relations_update, file="Gene_annotation/Gene_CpG_Relations_updatejune2015.RData")
