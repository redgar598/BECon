########## 
##BECon mQTL enrichment test
########## 
setwd("/big_data/redgar/BECon")




###### mQTL Data
ALSPEC_Middle_age<-read.table("FOM.ALL.M.tab", header=T)

#from http://www.mqtldb.org/download.htm
#the PLINK database reports associations (regardless of p-value) provided the SNP-CpG association was reported below 1e-7 in at least one timepoint;
## After conservative multiple testing correction (p<1×10−14) we identified between 24,262 and 31,729 sentinel associations at each time point

ALSPEC_Middle_age_sig<-ALSPEC_Middle_age[which(ALSPEC_Middle_age$p.value<=1e-14),]
#how many CpGs
mQTL_CpGs<-as.character(unique(ALSPEC_Middle_age_sig$gene))
length(mQTL_CpGs)




### Informative CpG analysis
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

neg<-unique(unlist(negCor_CpGs))#16738
pos<-unique(unlist(Cor_CpGs))#24453


informative<-unique(c(pos,neg))#40,029




########## 
## ENRICHMENT
########## 
informative_mQTL<-informative[which(informative%in%mQTL_CpGs)] #8202/40029

permutation_number=10000

random_CpG_mQTL<-sapply(1:permutation_number, function(x){
  print(x)
  set.seed(x)
  random_CpGs<-sample(correlations_BLBR_forblood$CpG, 40029)
  random_mQTL<-random_CpGs[which(random_CpGs%in%mQTL_CpGs)]
  length(random_mQTL)
})


#How many of the null means are bigger than the observed value? That proportion would be the
#p-value for the null. We add a 1 to the numerator and denominator to account for misestimation of
#the p-value (for more details see Phipson and Smyth, Permutation P-values should never be zero).
#the proportion of permutations with larger difference

(sum(random_CpG_mQTL > length(informative_mQTL)) + 1) / (permutation_number + 1)


## Just positive
pos_informative_mQTL<-pos[which(pos%in%mQTL_CpGs)] #6271/24453
permutation_number=10000
random_CpG_mQTL_pos<-sapply(1:permutation_number, function(x){
  print(x)
  set.seed(x)
  random_CpGs<-sample(correlations_BLBR_forblood$CpG, 24453)
  random_mQTL<-random_CpGs[which(random_CpGs%in%mQTL_CpGs)]
  length(random_mQTL)
})
(sum(random_CpG_mQTL_pos > length(pos_informative_mQTL)) + 1) / (permutation_number + 1)

## Just negative
neg_informative_mQTL<-neg[which(neg%in%mQTL_CpGs)] #2076/16738
permutation_number=10000
random_CpG_mQTL_neg<-sapply(1:permutation_number, function(x){
  print(x)
  set.seed(x)
  random_CpGs<-sample(correlations_BLBR_forblood$CpG, 16738)
  random_mQTL<-random_CpGs[which(random_CpGs%in%mQTL_CpGs)]
  length(random_mQTL)
})
(sum(random_CpG_mQTL_neg > length(neg_informative_mQTL)) + 1) / (permutation_number + 1)


########## 
## Hannon informative
##########
load("Positively_correlated_CpGs_JMill_atreference_range_0.1_corfrom0.2_nosex_nosnp.RData")
load("Negatively_correlated_CpGs_JMill_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")
load("GSE59685_correlations_blood_quantile_variation.RData")


neg<-unique(unlist(negCor_CpGs_JM))#42
pos<-unique(unlist(Cor_CpGs_JM))#10898

informative_JM<-unique(c(pos,neg))#10930


########## 
## ENRICHMENT
########## 
informative_mQTL_JM<-informative_JM[which(informative_JM%in%mQTL_CpGs)] #3018/10930

permutation_number=10000

random_CpG_mQTL<-sapply(1:permutation_number, function(x){
  print(x)
  set.seed(x)
  random_CpGs<-sample(correlations_BLBR_JMill$CpG, 10930)
  random_mQTL<-random_CpGs[which(random_CpGs%in%mQTL_CpGs)]
  length(random_mQTL)
})


#How many of the null means are bigger than the observed value? That proportion would be the
#p-value for the null. We add a 1 to the numerator and denominator to account for misestimation of
#the p-value (for more details see Phipson and Smyth, Permutation P-values should never be zero).
#the proportion of permutations with larger difference

(sum(random_CpG_mQTL > length(informative_mQTL_JM)) + 1) / (permutation_number + 1)


## Just positive
pos_informative_mQTL<-pos[which(pos%in%mQTL_CpGs)] #3013/10898
permutation_number=10000
random_CpG_mQTL_pos<-sapply(1:permutation_number, function(x){
  print(x)
  set.seed(x)
  random_CpGs<-sample(correlations_BLBR_forblood$CpG, 10898)
  random_mQTL<-random_CpGs[which(random_CpGs%in%mQTL_CpGs)]
  length(random_mQTL)
})
(sum(random_CpG_mQTL_pos > length(pos_informative_mQTL)) + 1) / (permutation_number + 1)

## Just negative
neg_informative_mQTL<-neg[which(neg%in%mQTL_CpGs)] #5/42
permutation_number=10000
random_CpG_mQTL_neg<-sapply(1:permutation_number, function(x){
  print(x)
  set.seed(x)
  random_CpGs<-sample(correlations_BLBR_forblood$CpG, 42)
  random_mQTL<-random_CpGs[which(random_CpGs%in%mQTL_CpGs)]
  length(random_mQTL)
})
(sum(random_CpG_mQTL_neg > length(neg_informative_mQTL)) + 1) / (permutation_number + 1)

