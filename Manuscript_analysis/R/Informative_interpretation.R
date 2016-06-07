### Informative CpG analysis
setwd("~/Documents/Blood_Brain/")
load("cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     64
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta
load("SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta_adjusted[which(!(rownames(combat_BLBR_Beta_adjusted)%in%rownames(SnpatCpG))),]

#Load correlation Level with blood varibility
load("Correlations_BLBR_blood_quantile_variation.RData") #already filtered for SNPs     
load("Price_annotation.RData")
Chrsex<-annotation$TargetID[which(annotation$CHR%in%c("X","Y"))]#11648
correlations_BLBR_forblood<-correlations_BLBR_forblood[which(!(correlations_BLBR_forblood$CpG%in%Chrsex)),]#413466

##Load top correlated CpGs
load("Positively_correlated_CpGs_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")
load("Negatively_correlated_CpGs_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")

neg<-unique(unlist(negCor_CpGs))
pos<-unique(unlist(Cor_CpGs))

combined<-c(pos,neg)#41,191
dup<-combined[which(duplicated(combined))]
dup_cor<-correlations_BLBR_forblood[which(correlations_BLBR_forblood$CpG%in%dup),]## 5 Sites are positive in one tissue and negative in another

informative<-unique(c(pos,neg))#40,029


####### HAVENT RUN ALL BELOW



########################################### informative genes FOR GO
load("Gene_CpG_Relations_updatejune2015.RData")

Gene_CpG_Relations_update<-Gene_CpG_Relations_update[!duplicated(Gene_CpG_Relations_update), ]
Genes_correlated_CpGs<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%informative),]
length(unique(Genes_correlated_CpGs$Probe_ID))#3587 CpGs
length(unique(Genes_correlated_CpGs$gene))#2028 Genes
Genes_correlated_CpGs<-Genes_correlated_CpGs[!duplicated(Genes_correlated_CpGs),]
Genes_correlated_CpGs<-Genes_correlated_CpGs[!duplicated(Genes_correlated_CpGs[,c(1,4)]),]#remove duplicate CpG to gene associations
informative_CpGs_all<-Genes_correlated_CpGs



#Suprising?
#suprising wiht the number of probes in the gene?
## Gene CpG number adjustment
Gene_CpG_Relations_update$gene<-as.character(Gene_CpG_Relations_update$gene)
Overrep<-as.data.frame(tapply(Gene_CpG_Relations_update$Probe_ID, Gene_CpG_Relations_update$gene, length))
Overrep$Gene<-rownames(Overrep)
colnames(Overrep)<-c("CpG_number", "Gene")
Overrep<-Overrep[which(Overrep$Gene!="None"),]
mean(Overrep$CpG_number, na.rm=T)# 25
Overrep$Enrichment_fromAverage<-Overrep$CpG_number/mean(Overrep$CpG_number, na.rm=T)

## Gene summaries
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

informative_CpGs_genes_all_multihits<-informative_CpGs_genes_all[which(informative_CpGs_genes_all$CpG_Associated>10),]
write.csv(as.character(unique(informative_CpGs_genes_all_multihits$gene)), file="Informative_CpGs_MultipleHitGenes_pos_neg_allchr.csv", row.names=FALSE)



##ORA results
ora<-read.csv("ORA_run.erminej_multiHits_informative10CPGS.csv", sep="\t")

######################################################## COEXPRESSION

## random gene lists, need to repeat with seeds
set.seed(1)
random_gene_list<-unique(Gene_CpG_Relations_update$gene)[sample(1:length(unique(Gene_CpG_Relations_update$gene)),188)]
write.csv(random_gene_list, file="Informative_CpGs_random_gene_list.csv")

set.seed(2)
random_gene_list<-unique(Gene_CpG_Relations_update$gene)[sample(1:length(unique(Gene_CpG_Relations_update$gene)),188)]
write.csv(random_gene_list, file="Informative_CpGs_random2_gene_list.csv")

set.seed(3)
random_gene_list<-unique(Gene_CpG_Relations_update$gene)[sample(1:length(unique(Gene_CpG_Relations_update$gene)),188)]
write.csv(random_gene_list, file="Informative_CpGs_random3_gene_list.csv")

set.seed(4)
random_gene_list<-unique(Gene_CpG_Relations_update$gene)[sample(1:length(unique(Gene_CpG_Relations_update$gene)),188)]
write.csv(random_gene_list, file="Informative_CpGs_random4_gene_list.csv")

set.seed(5)
random_gene_list<-unique(Gene_CpG_Relations_update$gene)[sample(1:length(unique(Gene_CpG_Relations_update$gene)),188)]
write.csv(random_gene_list, file="Informative_CpGs_random5_gene_list.csv")



########## co expression
co_expression_correlated<-read.csv("~/Documents/Blood_Brain/Gemma_coexpression/Informative_multipleHits_allChr", skip=9, sep="\t")
mean(co_expression_correlated$Positive.Support)
nrow(co_expression_correlated)
co_expression_correlated$Dataset<-"Informative CpGs"

co_expression_rnd1<-read.csv("~/Documents/Blood_Brain/Gemma_coexpression/Informative_random_allChr", skip=9, sep="\t")
mean(co_expression_rnd1$Positive.Support)
nrow(co_expression_rnd1)
co_expression_rnd1$Dataset<-"Random1"

co_expression_rnd2<-read.csv("~/Documents/Blood_Brain/Gemma_coexpression/Informative_random_allChr1", skip=9, sep="\t")
mean(co_expression_rnd2$Positive.Support)
nrow(co_expression_rnd2)
co_expression_rnd2$Dataset<-"Random2"

co_expression_rnd3<-read.csv("~/Documents/Blood_Brain/Gemma_coexpression/Informative_random_allChr2", skip=9, sep="\t")
mean(co_expression_rnd3$Positive.Support)
nrow(co_expression_rnd3)
co_expression_rnd3$Dataset<-"Random3"

co_expression_rnd4<-read.csv("~/Documents/Blood_Brain/Gemma_coexpression/Informative_random_allChr3", skip=9, sep="\t")
mean(co_expression_rnd4$Positive.Support)
nrow(co_expression_rnd4)
co_expression_rnd4$Dataset<-"Random4"

co_expression_rnd5<-read.csv("~/Documents/Blood_Brain/Gemma_coexpression/Informative_random_allChr4", skip=9, sep="\t")
mean(co_expression_rnd5$Positive.Support)
nrow(co_expression_rnd5)
co_expression_rnd5$Dataset<-"Random5"



coexp<-rbind(co_expression_correlated,co_expression_rnd1,co_expression_rnd2,co_expression_rnd3,co_expression_rnd4,co_expression_rnd5)

tapply(coexp$Positive.Support,coexp$Dataset, mean)
(tapply(coexp$Positive.Support,coexp$Dataset, mean)/tapply(coexp$Datasets.tested,coexp$Dataset, mean))*100
tapply(coexp$Negative.Support,coexp$Dataset, mean)
tapply(coexp$Query.Gene,coexp$Dataset, length)


coexp$percent_pos<-(coexp$Positive.Support/coexp$Datasets.tested)*100
ggplot(coexp, aes(Dataset, percent_pos))+geom_boxplot()

#informative CpGs have 179 connections vs mn of 55 connection in random lists of 188 genes.
# also informative CpGs assocaite with gene with slightly higher positive support for coexpression





################################################ Protein Protein interaction score
ReactomeF1<-read.csv("~/Documents/Blood_Brain/Reactome_protein_protein/FIsInGene_121013_with_annotations.txt", sep="\t")
Genes_correlated_multihits<-unique(Genes_correlated_CpGs_genes_all_multihits$gene)#120

## PP interactions between correlated CpGs
correlated_PP<-ReactomeF1[which(ReactomeF1$Gene1%in%Genes_correlated_multihits),]
correlated_PP<-correlated_PP[which(correlated_PP$Gene2%in%Genes_correlated_multihits),]#22 PP interactions

random_PP_number<-sapply(1:10000, function(x){
  set.seed(x)
  random_gene_list<-unique(Gene_CpG_Relations_update$gene)[sample(1:length(unique(Gene_CpG_Relations_update$gene)),188)]
  Genes_rnd_PP<-ReactomeF1[which(ReactomeF1$Gene1%in%random_gene_list),]
  Genes_rnd_PP<-Genes_rnd_PP[which(Genes_rnd_PP$Gene2%in%random_gene_list),]#4 PP interactions
  nrow(Genes_rnd_PP)})

ggplot()+geom_density(aes(random_PP_number), fill="grey")+geom_vline(xintercept=nrow(correlated_PP))+theme_bw()



############################################### Are informative sites closer to eachother than expected
annotation$CHR<-as.character(annotation$CHR)

#split by chr
all_annotation<-split(annotation, annotation$CHR)

each_Chr_all<-lapply(1:length(all_annotation), function(chr){
  annotation_df<-all_annotation[[chr]]
  distance_to_closest_CpG<-sapply(1:nrow(annotation_df), function(x){
    coor<-annotation_df$MAPINFO[x]
    coor_closest<-annotation_df$MAPINFO[-x][which.min(abs(annotation_df$MAPINFO[-x] - coor)) ]
    coor-coor_closest})})

invisible(lapply(1:length(all_annotation), function(chr) all_annotation[[chr]]$distance_to_closest_CpG<<-each_Chr_all[[chr]]))
invisible(lapply(1:length(all_annotation), function(chr) all_annotation[[chr]]$CpG<<-rownames(all_annotation[[chr]])))

all_annotation<-do.call(rbind, all_annotation)

save(all_annotation, file="Price_annotation_with_distance_to_closest_CpG.RData")

# parse by pos and neg from informative
pos_annotation<-all_annotation[which(all_annotation$CpG%in%pos),]
neg_annotation<-all_annotation[which(all_annotation$CpG%in%neg),]
informative_annotation<-all_annotation[which(all_annotation$CpG%in%combined),]


summary(abs(pos_annotation$distance_to_closest_CpG))
summary(abs(neg_annotation$distance_to_closest_CpG))
summary(abs(informative_annotation$distance_to_closest_CpG))
summary(abs(all_annotation$distance_to_closest_CpG))


pos_annotation$direction<-"Positive"
neg_annotation$direction<-"Negative"
all_annotation$direction<-"All"

direction_plot<-rbind(pos_annotation, neg_annotation, all_annotation)

ggplot(direction_plot[which(abs(direction_plot$distance_to_closest_CpG)<200000),], aes(direction, log(abs(distance_to_closest_CpG))))+
  geom_violin(fill="lightgrey", color="white")+geom_boxplot(width=0.25, fill="lightblue", outlier.size=0.5, outlier.shape=19)+theme_bw()

ggplot(direction_plot[which(abs(direction_plot$distance_to_closest_CpG)<200000),], aes(direction, (abs(distance_to_closest_CpG))))+
  geom_violin(fill="lightgrey", color="white")+geom_boxplot(width=0.25, fill="lightblue", outlier.size=0.5, outlier.shape=19)+theme_bw()


wilcox.test(log(pos_annotation$distance_to_closest_CpG), log(neg_annotation$distance_to_closest_CpG))
wilcox.test(pos_annotation$distance_to_closest_CpG, all_annotation$distance_to_closest_CpG)
wilcox.test(neg_annotation$distance_to_closest_CpG, all_annotation$distance_to_closest_CpG)





