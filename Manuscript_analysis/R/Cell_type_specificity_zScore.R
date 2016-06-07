### Celltype specificy z score

setwd("/big_data/redgar/echelon_backup_2/Blood_Brain")
library(reshape)
library(ggplot2)
library(RColorBrewer)


# After cell adjustment
load("cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     6
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta

#before adjustment
load("BLBR_Beta_Combatted_Mval_bmiqonly.RData") #441198     70
combat_BLBR_Beta_unadjusted<-as.data.frame(combat_BLBR_Beta)


# get them in the same order and porbe/sample number
combat_BLBR_Beta_unadjusted<-combat_BLBR_Beta_unadjusted[,which(colnames(combat_BLBR_Beta_unadjusted)%in%colnames(combat_BLBR_Beta_adjusted))]# 441198 n=63
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta_adjusted[,which(colnames(combat_BLBR_Beta_adjusted)%in%colnames(combat_BLBR_Beta_unadjusted))]
combat_BLBR_Beta_unadjusted<-combat_BLBR_Beta_unadjusted[which(rownames(combat_BLBR_Beta_unadjusted)%in%rownames(combat_BLBR_Beta_adjusted)),]# 441198 n=6
identical(rownames(combat_BLBR_Beta_unadjusted), rownames(combat_BLBR_Beta_adjusted))
identical(colnames(combat_BLBR_Beta_unadjusted), colnames(combat_BLBR_Beta_adjusted))

# Mean Delta Beta for blood before and after
meta<-read.csv("SampleInfo.csv")
fix_names<-unlist(lapply(meta$X, function(x) gsub(" ","", x , fixed=TRUE)))
meta$X<-fix_names
meta<-meta[which(meta$X%in%colnames(combat_BLBR_Beta_adjusted)),]
meta<-meta[match(colnames(combat_BLBR_Beta_adjusted), meta$X),]




####### Blood

# just blood 
PBMC_indx<-which(meta$TissueType=="PBMC") #16
combat_BL_Beta<-combat_BLBR_Beta_unadjusted[,PBMC_indx]
combat_BL_Beta_adjusted<-combat_BLBR_Beta_adjusted[,PBMC_indx]

Blood_changes<-do.call(rbind,lapply(1:nrow(combat_BL_Beta_adjusted), function(x) combat_BL_Beta[x,]-combat_BL_Beta_adjusted[x,]))
Blood_changes_abs_mn<-sapply(1:nrow(Blood_changes), function(x) mean(as.numeric((abs(Blood_changes[x,]))), na.rm=T))

mn_bld<-mean(Blood_changes_abs_mn)
sd_bld<-sd(Blood_changes_abs_mn)

Blood_changes_zscore<-sapply(1:length(Blood_changes_abs_mn), function(x) (Blood_changes_abs_mn[x]-mn_bld)/sd_bld)




####### Brain

# just brain 
brain_indx<-grep("BRAIN",meta$TissueType) #16
combat_BR_Beta<-combat_BLBR_Beta_unadjusted[,brain_indx]
combat_BR_Beta_adjusted<-combat_BLBR_Beta_adjusted[,brain_indx]

brain_changes<-do.call(rbind,lapply(1:nrow(combat_BR_Beta_adjusted), function(x) combat_BR_Beta[x,]-combat_BR_Beta_adjusted[x,]))
brain_changes_abs_mn<-sapply(1:nrow(brain_changes), function(x) mean(as.numeric((abs(brain_changes[x,]))), na.rm=T))

mn_brn<-mean(brain_changes_abs_mn)
sd_brn<-sd(brain_changes_abs_mn)

brain_changes_zscore<-sapply(1:length(brain_changes_abs_mn), function(x) (brain_changes_abs_mn[x]-mn_brn)/sd_brn)



## save
Brain_Blood_changes<-data.frame(CpG=rownames(combat_BL_Beta_adjusted),Blood_mean_beta_change=Blood_changes_abs_mn, Blood_z=Blood_changes_zscore, Brain_mean_beta_change=brain_changes_abs_mn, Brain_z=brain_changes_zscore)
save(Brain_Blood_changes, file="Brain_Blood_changes_zscore.Rdata")


load("Brain_Blood_changes_zscore.Rdata")
## want a beta value to display

## percentiles of delta beta change
Blood_percentile<-quantile(Brain_Blood_changes$Blood_mean_beta_change, c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
Brain_percentile<-quantile(Brain_Blood_changes$Brain_mean_beta_change, c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))

percential_blood<-sapply(Brain_Blood_changes$Blood_mean_beta_change, function(dB){
    if(dB<=Blood_percentile[5]){"<50% (Cell Composition Blood Percentile)"}else{
      if(dB<=Blood_percentile[7]){"50-75% (Cell Composition Blood Percentile)"}else{
        if(dB<=Blood_percentile[9]){"75-90% (Cell Composition Blood Percentile)"}else{"90% (Cell Composition Blood Percentile)"}}}}
)

percential_brain<-sapply(Brain_Blood_changes$Brain_mean_beta_change, function(dB){
  if(dB<=Brain_percentile[5]){"<50% (Cell Composition Brain Percentile)"}else{
    if(dB<=Brain_percentile[7]){"50-75% (Cell Composition Brain Percentile)"}else{
      if(dB<=Brain_percentile[9]){"75-90% (Cell Composition Brain Percentile)"}else{"90% (Cell Composition Brain Percentile)"}}}}
)

Brain_Blood_changes$percential_blood<-percential_blood
Brain_Blood_changes$percential_brain<-percential_brain

Brain_Blood_changes[,2]<-sapply(Brain_Blood_changes[,2], function(x) round(x, digits=2))
Brain_Blood_changes[,4]<-sapply(Brain_Blood_changes[,4], function(x) round(x, digits=2))



Brain_Blood_changes_percentile<-melt(Brain_Blood_changes[,c(1,6:7)], id="CpG")
Brain_Blood_changes_dB<-melt(Brain_Blood_changes[,c(1,2,4)], id=c("CpG"))

Brain_Blood_changes_melt<-cbind(Brain_Blood_changes_dB, Brain_Blood_changes_percentile)
Brain_Blood_changes_melt$Data<-"Cell_composition"
colnames(Brain_Blood_changes_melt)[6]<-"color"
Brain_Blood_changes_melt<-Brain_Blood_changes_melt[,c(1,7,2,3,6)]

save(Brain_Blood_changes_melt, file="Brain_Blood_changes_for_app.RData")

load("BLBR_app_Objects.RData")
CpG_gene_melt_combined<-rbind(CpG_gene_melt_combined, Brain_Blood_changes_melt)
