### Informative CpG analysis
setwd("/big_data/redgar/BECon")
load("cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     64
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta
load("SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta_adjusted[which(!(rownames(combat_BLBR_Beta_adjusted)%in%rownames(SnpatCpG))),]

#Load correlation levels and blood varibility metrics
load("Correlations_BLBR_blood_quantile_variation.RData")

# load 450K annotation
load("Price_annotation.RData")
Chrsex<-annotation$TargetID[which(annotation$CHR%in%c("X","Y"))]#11648
correlations_BLBR_forblood<-correlations_BLBR_forblood[which(!(correlations_BLBR_forblood$CpG%in%Chrsex)),]#413466

##Load informative CpGs both positivly and negatively correalated
load("Positively_correlated_CpGs_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")
load("Negatively_correlated_CpGs_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")

neg<-unique(unlist(negCor_CpGs))
pos<-unique(unlist(Cor_CpGs))

combined<-c(pos,neg)#41,191
dup<-combined[which(duplicated(combined))]
dup_cor<-correlations_BLBR_forblood[which(correlations_BLBR_forblood$CpG%in%dup),]## 5 Sites are positive in one tissue and negative in another

informative<-unique(c(pos,neg))#40,029


####### HAVENT RUN ALL BELOW


#############
## informative genes for GO enrichment analysis
#############
load("Gene_annotation/Gene_CpG_Relations_updatejune2015.RData")

Gene_CpG_Relations_update<-Gene_CpG_Relations_update[!duplicated(Gene_CpG_Relations_update), ]
Genes_correlated_CpGs<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%informative),]
length(unique(Genes_correlated_CpGs$Probe_ID))#3587 CpGs
length(unique(Genes_correlated_CpGs$gene))#2028 Genes
Genes_correlated_CpGs<-Genes_correlated_CpGs[!duplicated(Genes_correlated_CpGs),]
Genes_correlated_CpGs<-Genes_correlated_CpGs[!duplicated(Genes_correlated_CpGs[,c(1,4)]),]#remove duplicate CpG to gene associations
informative_CpGs_all<-Genes_correlated_CpGs



# Is it suprising to see that many CpGs in a gene given how many are on the 450k?
Gene_CpG_Relations_update$gene<-as.character(Gene_CpG_Relations_update$gene)
Overrep<-as.data.frame(tapply(Gene_CpG_Relations_update$Probe_ID, Gene_CpG_Relations_update$gene, length))
Overrep$Gene<-rownames(Overrep)
colnames(Overrep)<-c("CpG_number", "Gene")
Overrep<-Overrep[which(Overrep$Gene!="None"),]
mean(Overrep$CpG_number, na.rm=T)# 25
Overrep$Enrichment_fromAverage<-Overrep$CpG_number/mean(Overrep$CpG_number, na.rm=T)

## Gene summary table
load("Price_annotation.RData")
annotation$CpG<-rownames(annotation)
format_geneTable<-function(Genes_correlated_CpGs){
  print(paste("CpGs Associated: ", length(unique(Genes_correlated_CpGs$Probe_ID)), sep=""))
  print(paste("Genes Associated: ", length(unique(Genes_correlated_CpGs$gene)), sep=""))
  Overrep_subset<-as.data.frame(tapply(Genes_correlated_CpGs$Probe_ID, as.character(Genes_correlated_CpGs$gene), length))
  Overrep_subset$Gene<-rownames(Overrep_subset)
  colnames(Overrep_subset)<-c("CpG_number", "Gene")
  Overrep_subset<-Overrep_subset[which(Overrep_subset$Gene!="None"),]
  Overrep_subset_merge<-merge(Overrep_subset, Overrep, by="Gene")
  colnames(Overrep_subset_merge)<-c("Gene","CpG_Associated","CpG_in_Gene", "Enrichment_fromAverage")
  Overrep_subset_merge$Suprise<-Overrep_subset_merge$CpG_Associated/Overrep_subset_merge$Enrichment_fromAverage
  Gene_table<-merge(Genes_correlated_CpGs, Overrep_subset_merge, by.x="gene", by.y="Gene")
  Gene_table<-merge(Gene_table, annotation[,c(49,50,58)], by.x="Probe_ID", by.y="CpG")
  
  pval<-data.frame(CpG=rownames(combat_BLBR_Beta_adjusted)[which(rownames(combat_BLBR_Beta_adjusted)%in%Gene_table$Probe_ID)],
                   BA7_cor=correlations_BLBR_forblood$BRAIN7[which(rownames(combat_BLBR_Beta_adjusted)%in%Gene_table$Probe_ID)],
                   BA10_cor=correlations_BLBR_forblood$BRAIN10[which(rownames(combat_BLBR_Beta_adjusted)%in%Gene_table$Probe_ID)],
                   BA20_cor=correlations_BLBR_forblood$BRAIN20[which(rownames(combat_BLBR_Beta_adjusted)%in%Gene_table$Probe_ID)],
                   CV_blood=correlations_BLBR_forblood$CV_blood[which(rownames(combat_BLBR_Beta_adjusted)%in%Gene_table$Probe_ID)])
  
  Gene_table<-merge(Gene_table, pval, by.x="Probe_ID", by.y="CpG")
  Gene_table<-Gene_table[,c(2,7,8,9,10,1,4,3,6,5,11,12,13,14,15,16)]
  Gene_table<-Gene_table[order(-Gene_table$Suprise, Gene_table$gene),]
  Gene_table}

informative_CpGs_genes_all<-format_geneTable(informative_CpGs_all)
write.csv(Genes_correlated_CpGs_genes_all, file="Informative_CpGs_pos_neg_allchr.csv")

# highly informative (have 10 informative CpGs)
informative_CpGs_genes_all_multihits<-informative_CpGs_genes_all[which(informative_CpGs_genes_all$CpG_Associated>10),]
write.csv(as.character(unique(informative_CpGs_genes_all_multihits$gene)), file="Informative_CpGs_MultipleHitGenes_pos_neg_allchr.csv", row.names=FALSE)


##ORA results from Erminej
ora<-read.csv("ORA_run.erminej_multiHits_informative10CPGS.csv", sep="\t")

   
