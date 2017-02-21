setwd("/big_data/redgar/BECon/")
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(GEOmetadb)
library(sva)
library(limma)
library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kmanifest)
library(quadprog)
library(cets)
library(wateRmelon)
library(lumi)
library(mixtools)



##################
## Beta Data
##################
# Downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE59685&format=file&file=GSE59685%5Fbetas%2Ecsv%2Egz
GSE59685<-read.csv("GSE59685_betas.csv", skip=5)
GSE59685<-GSE59685[-1,]
x<-as.matrix(GSE59685)
y<-apply(x, 2, as.numeric)
GSE59685_numeric<-as.data.frame(y)
GSE59685_numeric<-GSE59685_numeric[,-1]
rownames(GSE59685_numeric)<-GSE59685$X
save(GSE59685_numeric, file="GSE59685_numeric.RData")



##################
##GEO meta data
##################library(GEOmetadb)
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
dbListFields(con, "gsm")

x<-dbGetQuery(con, "select title,series_id,gsm,source_name_ch1,characteristics_ch1 from gsm where gpl='GPL13534'")
meta<-x[which(x$series_id=="GSE59685"),]

meta$Subject<-unlist(sapply(1:nrow(meta), function(x){
  gsub("subjectid: ","",strsplit(meta$characteristics_ch1[x], ";\t")[[1]][1], fixed = T) # always field 2
}))
meta$barcode<-unlist(sapply(1:nrow(meta), function(x){
  gsub("barcode: ","",strsplit(meta$characteristics_ch1[x], ";\t")[[1]][2], fixed = T) # always field 2
}))
meta$ad.disease.status<-unlist(sapply(1:nrow(meta), function(x){
  gsub("ad.disease.status: ","",strsplit(meta$characteristics_ch1[x], ";\t")[[1]][3], fixed = T) # always field 2
}))
meta$braak.stage<-unlist(sapply(1:nrow(meta), function(x){
  gsub("braak.stage: ","",strsplit(meta$characteristics_ch1[x], ";\t")[[1]][4], fixed = T) # always field 2
}))
meta$Sex<-unlist(sapply(1:nrow(meta), function(x){
  gsub("Sex: ","",strsplit(meta$characteristics_ch1[x], ";\t")[[1]][5], fixed = T) # always field 2
}))
meta$age.blood<-unlist(sapply(1:nrow(meta), function(x){
  gsub("age.blood: ","",strsplit(meta$characteristics_ch1[x], ";\t")[[1]][6], fixed = T) # always field 2
}))
meta$age.brain<-unlist(sapply(1:nrow(meta), function(x){
  gsub("age.brain: ","",strsplit(meta$characteristics_ch1[x], ";\t")[[1]][7], fixed = T) # always field 2
}))
meta$Tissue<-unlist(sapply(1:nrow(meta), function(x){
  gsub("source tissue: ","",strsplit(meta$characteristics_ch1[x], ";\t")[[1]][8], fixed = T) # always field 2
}))

meta$description<-NULL
meta$source_name_ch1<-NULL
meta$characteristics_ch1<-NULL
meta$title<-NULL

Blood<-subset(meta, Tissue=="whole blood")
Brain<-subset(meta, Tissue!="whole blood")

Brain_matched<-Brain[which(Brain$Subject%in%Blood$Subject),]
Blood_matched<-Blood[which(Blood$Subject%in%Brain_matched$Subject),]
Meta_matched<-rbind(Brain_matched,Blood_matched)
tapply(Brain_matched$Subject, Brain_matched$Tissue, function(x) length(unique(x)))#71 individuals will all 4 brain regions

save(Meta_matched, file='Meta_matched_GSE59685.Rdata')

load("Meta_matched_GSE59685.Rdata")
GSE59685_matched<-GSE59685_numeric[,which(colnames(GSE59685_numeric)%in%Meta_matched$gsm)]
save(GSE59685_matched, file="GSE59685.RData")



##################
## Probe Filter
##################
load("GSE59685.RData")
load("Price_annotation.RData")
load("Meta_matched_GSE59685.Rdata")

Meta_matched<-Meta_matched[which(Meta_matched$gsm%in%colnames(GSE59685_matched)),]
Meta_matched<-Meta_matched[match(colnames(GSE59685_matched), Meta_matched$gsm),]



##### Remove Probes removed in original data
load("cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     64
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta
GSE59685_matched_filtered <- GSE59685_matched[rownames(GSE59685_matched)%in%rownames(combat_BLBR_Beta_adjusted), ] #  441198    366
annotation_filtered <- annotation[rownames(annotation)%in%rownames(GSE59685_matched_filtered), ]

##################
## Normalize
##################
    #The Csv file downloaded is dasen normalized
       

# Heat Plot Scree
source("BECon_scripts/Heat_scree_plot.R")
PCA_full<-princomp(GSE59685_matched_filtered[complete.cases(GSE59685_matched_filtered),]) # only using the 417619 complete rows for PCA
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)

Meta_matched$Sentrix_ID<-sapply(1:nrow(Meta_matched), function(x) strsplit(Meta_matched$barcode[x],"_")[[1]][1])
Meta_matched$Sentrix_pos<-sapply(1:nrow(Meta_matched), function(x) strsplit(Meta_matched$barcode[x],"_")[[1]][2])
Meta_matched$Sentrix_row<-sapply(1:nrow(Meta_matched), function(x) paste(strsplit(Meta_matched$Sentrix_pos[x],"")[[1]][2:3], sep="", collapse=""))

meta_categorical <- Meta_matched[, c(5,6,7,10,11,12,13)]  # input column numbers in meta that contain categorical variables
meta_continuous <- Meta_matched[, c(8,9)]  # input column numbers in meta that contain continuous variables

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

# how far do you want the plot to go?
PCs_to_view<-10

heat_scree_plot(Loadings, Importance)

##################
## Combat
##################
Mval<-function(beta) log2(beta/(1-beta))

## one sample removed as only sample on sentrix 7512560157
GSE59685_matched_filtered<-GSE59685_matched_filtered[,which(Meta_matched$Sentrix_ID!="7512560157")]# 365 samples
Meta_matched<-Meta_matched[which(Meta_matched$Sentrix_ID!="7512560157"),]

# GSE59685 data
pheno = Meta_matched
edata = apply(GSE59685_matched_filtered, 1, Mval) # need mvalues for combat
edata = as.data.frame(edata)
edata = t(edata)
pheno$Sentrix_ID<-as.factor(pheno$Sentrix_ID)

# Can not protect tissue because too confounded
#mod = model.matrix(~as.factor(Tissue), data=pheno)
#Sentrix_ID
batch = pheno$Sentrix_ID
combat_GSE59685_Beta = ComBat(dat=edata, batch=batch, mod=NULL, numCovs=NULL, par.prior=TRUE)
#Back to betas
betas<-function(M) 2^M/((2^M)+1)
combat_GSE59685_Beta = apply(combat_GSE59685_Beta, 1, betas) # need mvalues for combat
combat_GSE59685_Beta = as.data.frame(combat_GSE59685_Beta)
combat_GSE59685_Beta = t(combat_GSE59685_Beta)
combat_GSE59685_Beta<-combat_GSE59685_Beta
save(combat_GSE59685_Beta, file="GSE59685_filtered_combat.RData")

#heat scree again
PCA_full<-princomp(combat_GSE59685_Beta[complete.cases(combat_GSE59685_Beta),]) # only using the 417619 complete rows for PCA
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
meta_categorical <- Meta_matched[, c(5,6,7,10,11,12,13)]  # input column numbers in meta that contain categorical variables
meta_continuous <- Meta_matched[, c(8,9)]  # input column numbers in meta that contain continuous variables

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

heat_scree_plot(Loadings, Importance)

##################
## Cell Type Correct
##################

RGset = as.data.frame(combat_GSE59685_Beta) #  441198    365
source("BECon_scripts/ECC2.R")
blood_GSE<-which(Meta_matched$Tissue=="whole blood")
cellprop2 <- ECC2(RGset[,blood_GSE])
cellprop<-as.data.frame(cellprop2$counts)
cellprop$Sample<-rownames(cellprop)
save(cellprop, file="GSE59685_Bloodcellcomposition.RData")

# PLot
Blood_meta<-Meta_matched[which(Meta_matched$Tissue=="whole blood"),]
cellprop_melt<-melt(cellprop, id="Sample")
ggplot(cellprop_melt, aes(variable, value))+
  geom_point(shape=19, position = position_jitter(w = 0.1))+theme_bw()+
  xlab("Sample")+ylab("Cell Type Proportion")
## Megans code to correct for cell type
Blood_beta<-GSE59685_matched_filtered[,which(colnames(GSE59685_matched_filtered)%in%Blood_meta$gsm)]
# Impute NAs
  imputeMedianv3<-function(x) apply(x, 1, ,function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
Blood_beta_imputed<-t(imputeMedianv3(Blood_beta))
avebeta.lm<-apply(Blood_beta_imputed, 1, function(x){
  cellprop[colnames(Blood_beta_imputed),]->blood
  lm(x~CD8T+CD4T+NK+Bcell+Mono,data=blood)})
residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
colnames(residuals)<-colnames(Blood_beta_imputed)
adj.residuals<-residuals+matrix(apply(Blood_beta_imputed, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
save(adj.residuals, file="GSE59685_blood_celltype_adjusted.RData")

## Check how adjusted for cell type
RGset = as.data.frame(adj.residuals)
cellprop2 <- ECC2(RGset)
cellprop_adjusted<-as.data.frame(cellprop2$counts)
cellprop_adjusted$Sample<-rownames(cellprop_adjusted)
save(cellprop_adjusted, file="GSE59685_Bloodcellcomposition_afteradjustment.RData")



#### CETS

GSE59685_brain_betas<-GSE59685_matched_filtered[,which(Meta_matched$Tissue!="whole blood")]# 441198    291
meta_brain<-Meta_matched[which(Meta_matched$gsm%in%colnames(GSE59685_brain_betas)),]
meta_brain<-meta_brain[match(colnames(GSE59685_brain_betas), meta_brain$gsm),]
#Load calibration data set
# load "brain dataset" from data file in cetsBrain
load("BECon_scripts/cetsBrain/data/cetsBrain.rda") # click on cetsBrain.rda file to place in workspace
#Create the neuron and glia reference profiles:
modelIdx <- list(neuron = pdBrain$celltype == "N", glia = pdBrain$celltype ==  "G")
# getReference returns a 2-column matrix, representing reference profiles for the two cell types.
refProfile <- getReference(brain, modelIdx)
#Estimate the neuronal proportion:
prop <- estProportion(GSE59685_brain_betas, profile = refProfile)
prop<-as.data.frame(prop)
prop$glia<-apply(prop,1,function(x) 1-x)
colnames(prop)<- c("neuron", "glia")
save(prop, file="GSE59685_Braincellcomposition.RData")

# Plot
prop$sample<-rownames(prop)
porp_melt<-melt(prop)
prop_merge<-merge(porp_melt, meta_brain, by.x="sample", by.y="gsm")
ggplot(prop_merge, aes(variable, value, color=Tissue))+
  geom_point(shape=19, position=position_jitter(width=0.15))+
  theme_bw()+
  scale_color_manual(values=c("#fb6a4a","#ef3b2c","#cb181d","#cb181d"))+
  facet_wrap(~Tissue)
# fit methylation data for each probe in the dataset by the neuronal proportion
# Impute NAs
imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
GSE59685_brain_betas_imputed<-t(imputeMedianv3(GSE59685_brain_betas))
brainFC_GSE<-which(Meta_matched$Tissue=="frontal cortex")
brainEC_GSE<-which(Meta_matched$Tissue=="entorhinal cortex")
brainSTG_GSE<-which(Meta_matched$Tissue=="superior temporal gyrus")
brainCB_GSE<-which(Meta_matched$Tissue=="cerebellum")
# Split by Brain Region
FC_betas_imputed<-GSE59685_brain_betas_imputed[,brainFC_GSE]
EC_betas_imputed<-GSE59685_brain_betas_imputed[,brainEC_GSE]
STG_betas_imputed<-GSE59685_brain_betas_imputed[,brainSTG_GSE]
CB_betas_imputed<-GSE59685_brain_betas_imputed[,brainCB_GSE]
#Adjust Betas for cell types
Adjust<-function(BR){
  avebeta.lm<-apply(BR, 1, function(x){ ## FULL Adjustment!
    brain.sub<-prop[colnames(BR),]
    lm(x~neuron,data=brain.sub)})
  residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
  colnames(residuals)<-colnames(BR)
  adj.residuals<-residuals+matrix(apply(BR, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
  adj.residuals}
FC_adjusted<-Adjust(FC_betas_imputed)
EC_adjusted<-Adjust(EC_betas_imputed)
STG_adjusted<-Adjust(STG_betas_imputed)
CB_adjusted<-Adjust(CB_betas_imputed)
adj.residuals<-cbind(FC_adjusted, EC_adjusted, STG_adjusted, CB_adjusted)
adj.residuals<-adj.residuals[,match(meta_brain$gsm, colnames(adj.residuals))]
save(adj.residuals, file="GSE59685_Braincellcomposition_afteradjustment.RData")

# do composition again (on only complete cases, since NAs introduced by adjustment)
prop_adjusted<- estProportion(adj.residuals, profile = refProfile)
prop_adjusted<-as.data.frame(prop_adjusted)
prop_adjusted$glia<-apply(prop_adjusted,1,function(x) 1-x)
colnames(prop_adjusted)<- c("neuron", "glia")
save(prop_adjusted, file="GSE59685_Braincellcomposition_after_adjustment.RData") 


#### PLot of check cell compositions worked
#Brain
load("GSE59685_Braincellcomposition.RData")
load("GSE59685_Braincellcomposition_after_adjustment.RData")
prop_adjusted$sample<-rownames(prop_adjusted)
cellprop_adjustedmelt<-melt(prop_adjusted, id="sample")
cellprop_adjustedmelt$Data<-"adjusted"
cellprop_melt<-melt(prop, id="sample")
cellprop_melt$Data<-"unadjusted"
cellprop_melt_brain<-rbind(cellprop_melt, cellprop_adjustedmelt)
cellprop_melt_braintype<-merge(cellprop_melt_brain, meta_brain, by.x="sample", by.y="gsm")
#Blood
load("GSE59685_Bloodcellcomposition.RData")
load("GSE59685_Bloodcellcomposition_afteradjustment.RData")
GSE59685_blood_betas<-GSE59685_matched_filtered[,which(Meta_matched$Tissue=="whole blood")]# 441198    75
meta_blood<-Meta_matched[which(Meta_matched$gsm%in%colnames(GSE59685_blood_betas)),]
meta_blood<-meta_blood[match(colnames(GSE59685_blood_betas), meta_blood$gsm),]
cellprop_melt<-melt(cellprop, id="Sample")
cellprop_melt$Data<-"unadjusted"
cellprop_melt_adjusted<-melt(cellprop_adjusted, id="Sample")
cellprop_melt_adjusted$Data<-"adjusted"
cellprop_melt<-rbind(cellprop_melt, cellprop_melt_adjusted)
cellprop_melt<-merge(cellprop_melt, meta_blood, by.x="Sample", by.y="gsm")
colnames(cellprop_melt)[1]<-"sample"
cell_prop_melt<-rbind(cellprop_melt, cellprop_melt_braintype)
cell_prop_melt$Data<-as.factor(cell_prop_melt$Data)
cell_prop_melt$Data <- factor(cell_prop_melt$Data, levels = c("unadjusted", "adjusted"))

ggplot(cell_prop_melt, aes(variable, value, color=as.factor(variable)))+
  geom_point(shape=19, position = position_jitter(w = 0.2))+theme_bw()+xlab("Sample")+ylab("Cell Type Proportion")+
  facet_grid(Data~Tissue, scales="free_x")+
  scale_color_manual(values=c("#a50026", "#bd0026", "#e31a1c", "#fc4e2a", "#fd8d3c", "#feb24c", "#4575b4", "#74add1"),guide=FALSE)

##################
## Combine Data
##################
load(file="GSE59685_blood_celltype_adjusted.RData")
blood_adjusted<-as.data.frame(adj.residuals)
load(file="GSE59685_Braincellcomposition_afteradjustment.RData")
brain_adjusted<-as.data.frame(adj.residuals)
GSE59685<-cbind(blood_adjusted, brain_adjusted)# 441198    366
save(GSE59685, file="GSE59685_cellTypeAdjusted.RData")


##################
## Load matched data
##################
load("~/BLBR/Meta_matched_GSE59685.Rdata")
load("~/BLBR/GSE59685_cellTypeAdjusted.RData")
Meta_matched<-Meta_matched[which(Meta_matched$gsm%in%colnames(GSE59685)),]
Meta_matched<-Meta_matched[match(colnames(GSE59685), Meta_matched$gsm),]# 366  10

Meta_matched$ID<-sapply(1:nrow(Meta_matched), function(x) if(Meta_matched$Tissue[x]=="frontal cortex"){"FC"}else{
  if(Meta_matched$Tissue[x]=="cerebellum"){"CB"}else{
    if(Meta_matched$Tissue[x]=="entorhinal cortex"){"EC"}else{
      if(Meta_matched$Tissue[x]=="superior temporal gyrus"){"STG"}else{if(Meta_matched$Tissue[x]=="whole blood"){"WB"}}
    }
  }
})

Meta_matched$ID<-sapply(1:nrow(Meta_matched), function(x) paste(Meta_matched$Subject[x], Meta_matched$ID[x], sep="_"))

#load BECon data
load("cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     64
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta

#### Any meta data
meta<-read.csv("SampleInfo.csv")
fix_names<-unlist(lapply(meta$X, function(x) gsub(" ","", x , fixed=TRUE)))
meta$X<-fix_names
meta<-meta[which(meta$X%in%colnames(combat_BLBR_Beta)),]
meta<-meta[match(colnames(combat_BLBR_Beta), meta$X),]


## impute missing samples
blood_subjects<-Meta_matched$Subject[which(Meta_matched$Tissue=="whole blood")]
Meta_matched<-Meta_matched[which(Meta_matched$Subject%in%blood_subjects),]
GSE59685<-GSE59685[,which(colnames(GSE59685)%in%Meta_matched$gsm)]
Meta_matched<-Meta_matched[match(colnames(GSE59685), Meta_matched$gsm),]# 362  10

#cerebellum
cerebellum_subjects<-Meta_matched$Subject[which(Meta_matched$Tissue=="cerebellum")]
missCB<-blood_subjects[which(!(blood_subjects%in%cerebellum_subjects))]#4

EC_subjects<-Meta_matched$Subject[which(Meta_matched$Tissue=="entorhinal cortex")]
missEC<-blood_subjects[which(!(blood_subjects%in%EC_subjects))]#3

FC_subjects<-Meta_matched$Subject[which(Meta_matched$Tissue=="frontal cortex")]
missFC<-blood_subjects[which(!(blood_subjects%in%FC_subjects))]#1

STG_subjects<-Meta_matched$Subject[which(Meta_matched$Tissue=="superior temporal gyrus")]
missSTG<-blood_subjects[which(!(blood_subjects%in%STG_subjects))]#0

missing_meta<-data.frame(series_id="GSE59685",gsm=c(paste("CB_",missCB, sep=""), paste("EC_",missEC, sep=""), paste("FC_",missFC, sep="")),
                         Subject=c(missCB, missEC, missFC), 
                         barcode=NA,ad.disease.status=NA, braak.stage=NA, Sex=NA, age.blood=NA,age.brain=NA,
                         Tissue=c(rep("cerebellum", length(missCB)),
                                  rep("entorhinal cortex", length(missEC)),
                                  rep("frontal cortex", length(missFC))),
                         ID=c(paste(missCB,"CB", sep="_"), paste(missEC,"EC", sep="_"), paste(missFC,"FC", sep="_")))

missing_beta<-data.frame(CpG=rownames(GSE59685), CB_91=NA,CB_59=NA, CB_1=NA,CB_5=NA, EC_113=NA, EC_15=NA, EC_106=NA, FC_106=NA)
missing_beta$CpG<-NULL

GSE59685_matched<-cbind(GSE59685, missing_beta)
Meta_matched<-rbind(Meta_matched, missing_meta)
Meta_matched<-Meta_matched[which(Meta_matched$gsm%in%colnames(GSE59685_matched)),]
Meta_matched<-Meta_matched[match(colnames(GSE59685_matched), Meta_matched$gsm),]# 370  10
colnames(GSE59685_matched)<-Meta_matched$ID

save(GSE59685_matched, Meta_matched, file="JMill_betas_metas_cell_corrected_combatted.RData")


##################
## Correlation calculations between brains and blood
##################
load("JMill_betas_metas_cell_corrected_combatted.RData")

#### prep GSE data for python
blood_indx<-which(Meta_matched$Tissue=="whole blood") 
cerebellum_indx<-which(Meta_matched$Tissue=="cerebellum") 
EC_indx<-which(Meta_matched$Tissue=="entorhinal cortex") 
FC_indx<-which(Meta_matched$Tissue=="frontal cortex") 
STG_indx<-which(Meta_matched$Tissue=="superior temporal gyrus") 

blood_BLBR<-GSE59685_matched[,blood_indx]
cerebellum_BLBR<-GSE59685_matched[,cerebellum_indx]
EC_BLBR<-GSE59685_matched[,EC_indx]
FC_BLBR<-GSE59685_matched[,FC_indx]
STG_BLBR<-GSE59685_matched[,STG_indx]

# put in blood sample order
blood<-Meta_matched$Subject[which(Meta_matched$Tissue=="whole blood")]
cb<-Meta_matched$Subject[which(Meta_matched$Tissue=="cerebellum")]
ec<-Meta_matched$Subject[which(Meta_matched$Tissue=="entorhinal cortex")]
fc<-Meta_matched$Subject[which(Meta_matched$Tissue=="frontal cortex")]
stg<-Meta_matched$Subject[which(Meta_matched$Tissue=="superior temporal gyrus")]

blood_BLBR<-blood_BLBR[,order(blood)]
cerebellum_BLBR<-cerebellum_BLBR[,order(cb)]
EC_BLBR<-EC_BLBR[,order(ec)]
FC_BLBR<-FC_BLBR[,order(fc)]
STG_BLBR<-STG_BLBR[,order(stg)]

# for python ran GSE59685_Correlation_Oct19.py
write.csv(Meta_matched, file="Python/Meta_matched_BLBR_Oct19.csv")
write.csv(blood_BLBR, file="Python/blood_BLBR_Oct19.csv")
write.csv(cerebellum_BLBR, file="Python/cerebellum_BLBR_Oct19.csv")
write.csv(EC_BLBR, file="Python/EC_BLBR_Oct19.csv")
write.csv(FC_BLBR, file="Python/FC_BLBR_Oct19.csv")
write.csv(STG_BLBR, file="Python/STG_BLBR_Oct19.csv")

#Load correlations from python
correlations_BLBR_JMill<-read.csv("Python/correlation_celladjusted_oct19.csv", sep="\t")
correlations_BLBR_JMill$X<-rownames(GSE59685_matched)
colnames(correlations_BLBR_JMill)<-c("CpG","CB","EC","FC","STG")
save(correlations_BLBR_JMill, file="GSE59685_correlations_celladjusted.RData")


##################
## Reference Range Varibility
##################
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}

## Blood
Varibility_JMILL_BLBR_Beta_blood<-sapply(1:nrow(blood_BLBR), function(y) Variation(as.numeric(blood_BLBR[y,])))# across all blood samples
correlations_BLBR_JMill$CV_blood<-Varibility_JMILL_BLBR_Beta_blood

save(correlations_BLBR_JMill, file="GSE59685_correlations_blood_quantile_variation.RData")


##################
## Sex and polymorphic CpG correlations
##################
## GSE59685 data
load("JMill_betas_metas_cell_corrected_combatted.RData")

# SNP CpGs
load("SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]### These are polymorphic CpGs!

# Sex Chr CpGs
load("Price_annotation.RData")
Chrsex<-annotation$TargetID[which(annotation$CHR%in%c("X","Y"))]#11648

GSE59685_matched_sex<-GSE59685_matched[which(rownames(GSE59685_matched)%in%Chrsex),]
GSE59685_matched_SNPs<-GSE59685_matched[which(rownames(GSE59685_matched)%in%rownames(SnpatCpG)),]
GSE59685_matched_sexSNPs<-rbind(GSE59685_matched_SNPs, GSE59685_matched_sex)
GSE59685_matched_sexSNPs<-GSE59685_matched_sexSNPs[!duplicated(rownames(GSE59685_matched_sexSNPs)),]#27862

#### prep GSE data for python
blood_indx<-which(Meta_matched$Tissue=="whole blood") 
cerebellum_indx<-which(Meta_matched$Tissue=="cerebellum") 
EC_indx<-which(Meta_matched$Tissue=="entorhinal cortex") 
FC_indx<-which(Meta_matched$Tissue=="frontal cortex") 
STG_indx<-which(Meta_matched$Tissue=="superior temporal gyrus") 

blood_BLBR<-GSE59685_matched_sexSNPs[,blood_indx]
cerebellum_BLBR<-GSE59685_matched_sexSNPs[,cerebellum_indx]
EC_BLBR<-GSE59685_matched_sexSNPs[,EC_indx]
FC_BLBR<-GSE59685_matched_sexSNPs[,FC_indx]
STG_BLBR<-GSE59685_matched_sexSNPs[,STG_indx]

# put in blood sample order
blood<-Meta_matched$Subject[which(Meta_matched$Tissue=="whole blood")]
cb<-Meta_matched$Subject[which(Meta_matched$Tissue=="cerebellum")]
ec<-Meta_matched$Subject[which(Meta_matched$Tissue=="entorhinal cortex")]
fc<-Meta_matched$Subject[which(Meta_matched$Tissue=="frontal cortex")]
stg<-Meta_matched$Subject[which(Meta_matched$Tissue=="superior temporal gyrus")]

blood_BLBR<-blood_BLBR[,order(blood)]
cerebellum_BLBR<-cerebellum_BLBR[,order(cb)]
EC_BLBR<-EC_BLBR[,order(ec)]
FC_BLBR<-FC_BLBR[,order(fc)]
STG_BLBR<-STG_BLBR[,order(stg)]

# for python ran GSE59685_Correlation_snpsex.py
write.csv(blood_BLBR, file="Python/blood_BLBR_snpsex_Oct19.csv")
write.csv(cerebellum_BLBR, file="Python/cerebellum_BLBR_snpsex_Oct19.csv")
write.csv(EC_BLBR, file="Python/EC_BLBR_snpsex_Oct19.csv")
write.csv(FC_BLBR, file="Python/FC_BLBR_snpsex_Oct19.csv")
write.csv(STG_BLBR, file="Python/STG_BLBR_snpsex_Oct19.csv")

#Load correlations from python
correlations_BLBR_JMill_snpsex<-read.csv("Python/correlation_snpsex_celladjusted_oct19.csv", sep="\t")
correlations_BLBR_JMill_snpsex$X<-rownames(GSE59685_matched_sexSNPs)
colnames(correlations_BLBR_JMill_snpsex)<-c("CpG","CB","EC","FC","STG")
save(correlations_BLBR_JMill_snpsex, file="GSE59685_correlations_snpsex_celladjusted.RData")



##################
## Polymorphic and Sex ref range
##################
load("GSE59685_correlations_snpsex_celladjusted.RData")
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}

## Blood
Varibility_snpsex_JMILL_BLBR_Beta_blood<-sapply(1:nrow(blood_BLBR), function(y) Variation(as.numeric(blood_BLBR[y,])))# across all blood samples
correlations_BLBR_JMill_snpsex$CV_blood<-Varibility_snpsex_JMILL_BLBR_Beta_blood
save(correlations_BLBR_JMill_snpsex, file="GSE59685_correlations_snpsex_blood_quantile_variation.RData")


##################
## Compare Correlation DIstributions
##################
load("GSE59685_correlations_blood_quantile_variation.RData")
load("Correlations_BLBR_blood_quantile_variation.RData") #already filtered for SNPs             

#filter J Mill for SNPs
correlations_BLBR_JMill<-correlations_BLBR_JMill[which(correlations_BLBR_JMill$CpG%in%correlations_BLBR_forblood$CpG),]


### plot
plotcorr<-cbind(correlations_BLBR_forblood[,1:4], correlations_BLBR_JMill[,2:5])
plotcorr<-melt(plotcorr)

levels(plotcorr$variable)<-c("BA7","BA10","BA20","CB","EC","FC","STG")

ggplot(plotcorr, aes(value, color=variable))+geom_density(size=1)+theme_bw()+
  scale_color_manual(values=c("#fb6a4a","#ef3b2c","#cb181d","#810f7c","#88419d","#8c6bb1","#8c96c6"), name="Brain \nRegion")+
  xlab("Correlation Coefficient")


##################
## Mixture model with Polymorphic Sex CpGs
##################
load("GSE59685_correlations_snpsex_blood_quantile_variation.RData")
correlations_BLBR_JMill_ALL<-rbind(correlations_BLBR_JMill_snpsex, correlations_BLBR_JMill)
Cutoffs<-c(0.05, 0.075,0.1,0.125,0.15,0.175,0.2)
DensityPLot_snp_sex<-lapply(Cutoffs, function(x) {
  varible<-correlations_BLBR_JMill_ALL[which(correlations_BLBR_JMill_ALL$CV_blood>x),]
  varible$CV_cutoff<-x
  varible})
N_rnd<-sapply(1:length(Cutoffs), function(x) nrow(DensityPLot_snp_sex[[x]]))
DensityPLot_snp_sex<-lapply(1:length(Cutoffs), function(x) melt(DensityPLot_snp_sex[[x]],id=c("CpG","CV_blood","CV_cutoff")))
DensityPLot_snp_sex<-do.call(rbind, DensityPLot_snp_sex)
# Plot
col<-colorRampPalette(brewer.pal(9,"Blues"))(length(Cutoffs))
ggplot(DensityPLot_snp_sex, aes(value, color=as.factor(CV_cutoff)))+geom_density()+
  theme_bw()+facet_wrap(~variable, ncol=1)+scale_color_manual(values=col, name="Variance Cutoff")+
  xlab("Correlation")+ylab("CpG Density")+ylim(0,3)

DensityPLot_snp_sex<-lapply(Cutoffs, function(x) {
  varible<-correlations_BLBR_JMill_ALL[which(correlations_BLBR_JMill_ALL$CV_blood>x),]
  varible$CV_cutoff<-x
  varible})
mixMod<-lapply(1:length(DensityPLot_snp_sex), function(cv) lapply(2:5, function(x) normalmixEM(DensityPLot_snp_sex[[cv]][,x],k=2, maxit=500)))  # takes awhile 
save(mixMod,file="mixMod_JMill_referenceReange_SNPs_Xchr.RData")

load("mixMod_JMill_referenceReange_SNPs_Xchr.RData")
cv<-7
s<-4
pos_cor_pop<-which(mixMod[[cv]][[s]]$mu==max(mixMod[[cv]][[s]]$mu))
plot(mixMod[[cv]][[s]],which=2,xlab2="Correlation")
lines(density(DensityPLot_blood[[cv]][,(s+1)]), lty=2, lwd=3)
abline(v= mixMod[[cv]][[s]]$mu[pos_cor_pop]-(2*mixMod[[cv]][[s]]$sigma[pos_cor_pop]), col="blue", lwd=3)

cor_threshold<-do.call(rbind,lapply(1:length(mixMod), function(cv) sapply(1:4, function(s) {
  pos_cor_pop<-which(mixMod[[cv]][[s]]$mu==max(mixMod[[cv]][[s]]$mu))#select the highly positively correalted sub population(ie max mu)
  mixMod[[cv]][[s]]$mu[pos_cor_pop]-(2*mixMod[[cv]][[s]]$sigma[pos_cor_pop])})))#what cor value is mu-2sd

cor_threshold<-as.data.frame(cor_threshold)
colnames(cor_threshold)<-c("CB","EC","FC","STG")
cor_threshold$CV_cutoff<-Cutoffs


# CpG Numbers (Cor at 0.2 and reference range at 0.1)
cor_threshold[7,]#thresholds
sapply(1:4, function(br) length(correlations_BLBR_JMill$CpG[which(correlations_BLBR_JMill[,(br+1)]>=cor_threshold[7,br])]))#CpGs Passing
sapply(1:4, function(br) length(correlations_BLBR_JMill$CpG[which(correlations_BLBR_JMill[,(br+1)]>=cor_threshold[7,br] & correlations_BLBR_JMill$CV_blood>=0.1)]))#CpGs Passing

#overlap in CpGs passing between three brain regions
Cor_CpGs_JM<-lapply(1:4, function(br) correlations_BLBR_JMill$CpG[which(correlations_BLBR_JMill[,(br+1)]>=cor_threshold[7,br] & correlations_BLBR_JMill$CV>=0.1)])#CpGs Passing
save(Cor_CpGs_JM, file="Positively_correlated_CpGs_JMill_atreference_range_0.1_corfrom0.2_nosex_nosnp.RData")

## negative
negCor_CpGs_JM<-lapply(1:4, function(br) correlations_BLBR_JMill$CpG[which(correlations_BLBR_JMill[,(br+1)]<=(-cor_threshold[7,br]) & correlations_BLBR_JMill$CV>=0.1)])#CpGs Passing
save(negCor_CpGs_JM, file="Negatively_correlated_CpGs_JMill_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")

##################         
## All Informative CpG Overlap
##################
load("Positively_correlated_CpGs_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")
load("Negatively_correlated_CpGs_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")
load("Negatively_correlated_CpGs_JMill_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")
load("Positively_correlated_CpGs_JMill_atreference_range_0.1_corfrom0.2_nosex_nosnp.RData")


overlap<-function(tissueA, tissueB){(length(intersect(tissueA, tissueB))/length(tissueA))*100}

## Positive
Cor_CpGs[4:7]<-Cor_CpGs_JM
names(Cor_CpGs)<-c("BA7","BA10","BA20","CB","EC","FC","STG")

overlapping_correlation<-do.call(rbind,lapply(1:length(Cor_CpGs), function(y){
  sapply(1:length(Cor_CpGs), function(x) overlap(Cor_CpGs[[y]], Cor_CpGs[[x]]))}))
colnames(overlapping_correlation)<-names(Cor_CpGs)
rownames(overlapping_correlation)<-names(Cor_CpGs)

### Is the overlap suprising?
Permutate_overlap<-function(probe_list1, probe_list2){
  len1<-length(probe_list1)
  len2<-length(probe_list2)
  rnd1<-correlations_BLBR_forblood$CpG[sample(1:nrow(correlations_BLBR_forblood), len1)]
  rnd2<-correlations_BLBR_forblood$CpG[sample(1:nrow(correlations_BLBR_forblood), len2)]
  (length(intersect(rnd1, rnd2))/length(rnd1))*100}

# Permutate
permutated_overlapping_correlation<-do.call(rbind,lapply(1:length(Cor_CpGs), function(y){
  sapply(1:length(Cor_CpGs), function(x) {
    expected<-sapply(1:100, function(seed){
      set.seed(seed)
      Permutate_overlap(Cor_CpGs[[y]],Cor_CpGs[[x]])})
    mean(expected)
  })
}))
colnames(permutated_overlapping_correlation)<-names(Cor_CpGs)
rownames(permutated_overlapping_correlation)<-names(Cor_CpGs)

plot_df<-as.data.frame(overlapping_correlation)
plot_df$TissueA<-names(Cor_CpGs)
plot_df<-melt(plot_df)

#add permutation expectation
plot_df_rnd<-as.data.frame(permutated_overlapping_correlation)
plot_df_rnd$TissueA<-names(Cor_CpGs)
plot_df_rnd<-melt(plot_df_rnd)
plot_df$expected<-plot_df_rnd$value

# Remove repeated comparisons (only want to see percents with respect to row)
plot_df<-plot_df[which(!(plot_df$TissueA==plot_df$variable)),]
plot_df<-plot_df[which(!(plot_df$TissueA=="BA7")),]
plot_df<-plot_df[which(!(plot_df$variable=="STG")),]

plot_df<-plot_df[-c(which(plot_df$variable=="BA20"& plot_df$TissueA=="BA10"),
                   which(plot_df$variable=="CB"& plot_df$TissueA%in%c("BA10","BA20")),
                   which(plot_df$variable=="EC"& plot_df$TissueA%in%c("BA10","BA20","CB")),
                   which(plot_df$variable=="FC"& plot_df$TissueA%in%c("BA10","BA20","CB","EC"))),]

plot_df$TissueA <- factor(plot_df$TissueA, levels = c("BA7","BA10","BA20","CB","EC","FC","STG"))
plot_df$variable <- factor(plot_df$variable, levels = c("BA7","BA10","BA20","CB","EC","FC","STG"))

levels(plot_df$TissueA) <- c("Brodmann\nBrain 7","Brodmann\nBrain 10","Brodmann\nBrain 20",
                             "Cerebellum","Entorhinal\nCortex","Frontal\nCortex","Superior\nTemporal\nGyrus")
levels(plot_df$variable) <- c("Brodmann\nBrain 7","Brodmann\nBrain 10","Brodmann\nBrain 20",
                             "Cerebellum","Entorhinal\nCortex","Frontal\nCortex","Superior\nTemporal\nGyrus")

ggplot(plot_df, aes(variable, TissueA, fill = value)) +
  geom_tile(color = "black",size=0.5) +
  geom_text(aes(label = paste(round(value,2),"\n","", sep="")), color="black")+theme_bw()+
  geom_text(aes(label = paste("","\n","(",round(expected, 2),")", sep="")), color="black", size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text = element_text(size=12))+
  scale_fill_continuous(low="#c6dbef", high="#4292c6", guide="none")+
  xlab("")+ylab("")




################# Negative (overlap is minimal and not shown in paper)
negCor_CpGs[4:7]<-negCor_CpGs_JM
names(negCor_CpGs)<-c("BA7","BA10","BA20","CB","EC","FC","STG")

sapply(1:length(negCor_CpGs), function(x) length(negCor_CpGs[[x]]))

overlapping_correlation_neg<-do.call(rbind,lapply(1:length(negCor_CpGs), function(y){
  sapply(1:length(negCor_CpGs), function(x) overlap(negCor_CpGs[[y]], negCor_CpGs[[x]]))}))
colnames(overlapping_correlation_neg)<-names(negCor_CpGs)
rownames(overlapping_correlation_neg)<-names(negCor_CpGs)


# Permutate
permutated_overlapping_correlation_neg<-do.call(rbind,lapply(1:length(negCor_CpGs), function(y){
  sapply(1:length(negCor_CpGs), function(x) {
    expected<-sapply(1:100, function(seed){
      set.seed(seed)
      Permutate_overlap(negCor_CpGs[[y]],negCor_CpGs[[x]])})
    mean(expected)
  })
}))
colnames(permutated_overlapping_correlation_neg)<-names(negCor_CpGs)
rownames(permutated_overlapping_correlation_neg)<-names(negCor_CpGs)


plot_df_neg<-as.data.frame(overlapping_correlation_neg)
plot_df_neg$TissueA<-names(negCor_CpGs)
plot_df_neg<-melt(plot_df_neg)

#add permutation expectation
plot_df_rnd_neg<-as.data.frame(permutated_overlapping_correlation_neg)
plot_df_rnd_neg$TissueA<-names(negCor_CpGs)
plot_df_rnd_neg<-melt(plot_df_rnd_neg)
plot_df_neg$expected<-plot_df_rnd_neg$value

plot_df_neg<-plot_df_neg[which(!(plot_df_neg$TissueA==plot_df_neg$variable)),]
plot_df_neg<-plot_df_neg[which(!(plot_df_neg$TissueA=="BA7")),]
plot_df_neg<-plot_df_neg[which(!(plot_df_neg$variable=="STG")),]

plot_df_neg<-plot_df_neg[-c(which(plot_df_neg$variable=="BA20"& plot_df_neg$TissueA=="BA10"),
                    which(plot_df_neg$variable=="CB"& plot_df_neg$TissueA%in%c("BA10","BA20")),
                    which(plot_df_neg$variable=="EC"& plot_df_neg$TissueA%in%c("BA10","BA20","CB")),
                    which(plot_df_neg$variable=="FC"& plot_df_neg$TissueA%in%c("BA10","BA20","CB","EC"))),]

plot_df_neg$TissueA <- factor(plot_df_neg$TissueA, levels = c("BA7","BA10","BA20","CB","EC","FC","STG"))
plot_df_neg$variable <- factor(plot_df_neg$variable, levels = c("BA7","BA10","BA20","CB","EC","FC","STG"))

levels(plot_df_neg$TissueA) <- c("Brodmann\nBrain 7","Brodmann\nBrain 10","Brodmann\nBrain 20",
                             "Cerebellum","Entorhinal\nCortex","Frontal\nCortex","Superior\nTemporal\nGyrus")
levels(plot_df_neg$variable) <- c("Brodmann\nBrain 7","Brodmann\nBrain 10","Brodmann\nBrain 20",
                              "Cerebellum","Entorhinal\nCortex","Frontal\nCortex","Superior\nTemporal\nGyrus")

ggplot(plot_df_neg, aes(variable, TissueA, fill = value)) +
  geom_tile(color = "black",size=0.5) +
  geom_text(aes(label = paste(round(value,2),"\n","", sep="")), color="black")+theme_bw()+
  geom_text(aes(label = paste("","\n","(",round(expected, 2),")", sep="")), color="black", size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text = element_text(size=12))+
  scale_fill_continuous(low="#c6dbef", high="#4292c6", guide="none")+
  xlab("")+ylab("")


  
  
  
  
  
      
################# mQTL PLot

load("cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     64
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta


load("Correlations_BLBR_blood_quantile_variation.RData") #already filtered for SNPs     
load("Price_annotation.RData")
Chrsex<-annotation$TargetID[which(annotation$CHR%in%c("X","Y"))]#11648
correlations_BLBR_forblood<-correlations_BLBR_forblood[which(!(correlations_BLBR_forblood$CpG%in%Chrsex)),]#413466

## informative sites BECon
load("Positively_correlated_CpGs_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")
load("Negatively_correlated_CpGs_atreference_range_0.1_corfrom0.25_nosex_nosnp.RData")

neg<-unique(unlist(negCor_CpGs))
pos<-unique(unlist(Cor_CpGs))
informative<-unique(c(pos,neg))#40,029


## pick top sites to plot
informative_correlation<-correlations_BLBR_forblood[which(correlations_BLBR_forblood$CpG%in%informative),]

## Pick CPgs to plot
x<-informative_correlation[which(abs(informative_correlation$BRAIN7)>0.8 & informative_correlation$CV_blood>0.45),]
x
# Cherry pick representative examples
#mqtls in theirs hidden in ours cg06417478
# Perfect cg25457927
# great cg27201625
# interesting cg02655711

CpGs<-c("cg06417478","cg25457927","cg27201625","cg23899408","cg04657146","cg08041448")

mqtl<-c("cg06417478","cg23899408","cg04657146")
notmqtl<-c("cg25457927","cg27201625","cg08041448")


### load Hannon beta values 
load("JMill_betas_metas_cell_corrected_combatted.RData")

GSE59685_mini<-GSE59685_matched[which(rownames(GSE59685_matched)%in%CpGs),]
blood_indx<-which(Meta_matched$Tissue=="whole blood") 
brain_indx<-which(Meta_matched$Tissue!="whole blood") 
GSE59685_mini_blood<-GSE59685_mini[,blood_indx]
GSE59685_mini_brain<-GSE59685_mini[,brain_indx]
GSE59685_mini_blood$CpG<-rownames(GSE59685_mini_blood)
GSE59685_mini_brain$CpG<-rownames(GSE59685_mini_brain)
GSE59685_mini_blood_melted<-melt(GSE59685_mini_blood)
GSE59685_mini_blood_melted<-merge(GSE59685_mini_blood_melted, Meta_matched, by.x="variable", by.y="ID")
colnames(GSE59685_mini_blood_melted)[3]<-"Blood_DNAm"
GSE59685_mini_brain_melted<-melt(GSE59685_mini_brain)
GSE59685_mini_brain_melted<-merge(GSE59685_mini_brain_melted, Meta_matched, by.x="variable", by.y="ID")
colnames(GSE59685_mini_brain_melted)[3]<-"Brain_DNAm"
GSE59685_mini_melted<-merge(GSE59685_mini_blood_melted, GSE59685_mini_brain_melted, by=c("Subject","CpG"))


## BECon data
meta<-read.csv("SampleInfo.csv")
fix_names<-unlist(lapply(meta$X, function(x) gsub(" ","", x , fixed=TRUE)))
meta$X<-fix_names
meta<-meta[which(meta$X%in%colnames(combat_BLBR_Beta_adjusted)),]
meta<-meta[match(colnames(combat_BLBR_Beta_adjusted), meta$X),]

BLBR_mini<-combat_BLBR_Beta_adjusted[which(rownames(combat_BLBR_Beta_adjusted)%in%CpGs),]
blood_indx<-which(meta$TissueType=="PBMC") 
brain_indx<-which(meta$TissueType!="PBMC") 
BLBR_mini_blood<-BLBR_mini[,blood_indx]
BLBR_mini_brain<-BLBR_mini[,brain_indx]
BLBR_mini_blood$CpG<-rownames(BLBR_mini_blood)
BLBR_mini_brain$CpG<-rownames(BLBR_mini_brain)
BLBR_mini_blood_melted<-melt(BLBR_mini_blood)
BLBR_mini_blood_melted<-merge(BLBR_mini_blood_melted, meta, by.x="variable", by.y="X")
colnames(BLBR_mini_blood_melted)[3]<-"Blood_DNAm"
BLBR_mini_brain_melted<-melt(BLBR_mini_brain)
BLBR_mini_brain_melted<-merge(BLBR_mini_brain_melted, meta, by.x="variable", by.y="X")
colnames(BLBR_mini_brain_melted)[3]<-"Brain_DNAm"

## Combined data for plotting
BLBR_mini_melted<-merge(BLBR_mini_blood_melted, BLBR_mini_brain_melted, by=c("SubjectNumC","CpG"))
BLBR_mini_melted<-BLBR_mini_melted[,c(1,2,4,17,22)]
GSE59685_mini_melted<-GSE59685_mini_melted[,c(1,2,4,15,24)]
colnames(BLBR_mini_melted)<-c("Subject","CpG","Blood_DNAm","Brain_DNAm","Tissue.y")
BLBR_mini_melted$Data<-"BLBR"
GSE59685_mini_melted$Data<-"GSE59685"
plot_data<-rbind(GSE59685_mini_melted, BLBR_mini_melted)

plot_data$Type<-sapply(1:nrow(plot_data), function(x) if(plot_data$CpG[x]%in%mqtl){"mQTL"}else{"Not mQTL"})

ggplot(plot_data, aes(Blood_DNAm, Brain_DNAm, fill=Tissue.y, color=Data, size=Data))+
  geom_point(shape=21)+theme_bw()+facet_wrap(Type~CpG)+scale_color_manual(values=c("black","white"), name="Dataset")+
  scale_fill_manual(values=c("#fb6a4a","#ef3b2c","#cb181d","#810f7c","#88419d","#8c6bb1","#8c96c6"), name="Brain \nRegion")+
  ylab("Brain Methylation")+xlab("Blood Methylation")+scale_size_manual(values=c(3,2))

  


