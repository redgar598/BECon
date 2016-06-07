#### validation genes from studies
setwd("~/Documents/Blood_Brain/")
load("Price_annotation.RData")
annotation$CpG<-rownames(annotation)
load("Gene_CpG_Relations_updatejune2015.RData")


## BDNF
BDNF<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$gene=="BDNF"),]
length(unique(BDNF$Probe_ID))


#A CpG-rich region of the promoter between −694 and −577, relative to the transcriptional start and including seven CpG sites
# has many variants which one did they use as TSS
BDNF_promoter<-BDNF[which(BDNF$region=="promoter"),]
length(unique(BDNF_promoter$Probe_ID))
#TSS chr11:27722600
# region 27722600
tss<-27722600
tss+557
tss+694

BDNF_promoter[which(BDNF_promoter$Coordinate_37>(tss+557) & BDNF_promoter$Coordinate_37<(tss+694)),]


###DRD2
chr11:113345999-113346540

DRD2<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$gene=="DRD2"),]
length(unique(DRD2$Probe_ID))

DRD2_promoter<-DRD2[which(DRD2$region=="promoter"),]
length(unique(DRD2_promoter$Probe_ID))
chr11<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Chromosome_37=="11"),]
chr11DRD2reg<-chr11[which(chr11$Coordinate_37<113346540 & chr11$Coordinate_37>113345999),]


###OXTR

chr3:8808689-8812491

OXTR<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$gene=="OXTR"),]
length(unique(OXTR$Probe_ID))

OXTR_promoter<-OXTR[which(OXTR$region=="promoter"),]
length(unique(OXTR_promoter$Probe_ID))
chr3<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Chromosome_37=="3"),]
chr3OXTRreg<-chr3[which(chr3$Coordinate_37<8812491 & chr3$Coordinate_37>8808689),]


###FAM63B
FAM63B<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$gene=="FAM63B"),]
length(unique(FAM63B$Probe_ID))

OXTR_promoter<-OXTR[which(OXTR$region=="promoter"),]
length(unique(OXTR_promoter$Probe_ID))
chr3<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Chromosome_37=="3"),]
chr3OXTRreg<-chr3[which(chr3$Coordinate_37<8812491 & chr3$Coordinate_37>8808689),]
