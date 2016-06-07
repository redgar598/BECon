setwd("~/Documents/Blood_Brain/")
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(Biobase)
library(mixtools)


#just the 16 for each sample and one brain20 NA
load("cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     64
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta
load("SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta_adjusted[which(!(rownames(combat_BLBR_Beta_adjusted)%in%rownames(SnpatCpG))),]

correlations_BLBR<-read.csv("~/Documents/Blood_Brain/Python/correlation_jan22.csv", sep="\t")                       
colnames(correlations_BLBR)<-c("CpG","BRAIN7","BRAIN10","BRAIN20")
correlations_BLBR$CpG<-rownames(combat_BLBR_Beta)
correlations_BLBR<-correlations_BLBR[which(!(correlations_BLBR$CpG%in%rownames(SnpatCpG))),]



########################## call varibility in each tissue
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}


## Just blood covariance
## Blood
blood<-grep("PBMC", colnames(combat_BLBR_Beta_adjusted))
allbrain<-grep("BA", colnames(combat_BLBR_Beta_adjusted))
BA7<-grep("BA7", colnames(combat_BLBR_Beta_adjusted))
BA10<-grep("BA10", colnames(combat_BLBR_Beta_adjusted))
BA20<-grep("BA20", colnames(combat_BLBR_Beta_adjusted))


Covariance_blood<-sapply(1:nrow(combat_BLBR_Beta_adjusted), function(y) Variation(as.numeric(combat_BLBR_Beta_adjusted[y,blood])))# across all blood samples
Covariance_allbrain<-sapply(1:nrow(combat_BLBR_Beta_adjusted), function(y) Variation(as.numeric(combat_BLBR_Beta_adjusted[y,allbrain])))# across all brain samples
Covariance_ba7<-sapply(1:nrow(combat_BLBR_Beta_adjusted), function(y) Variation(as.numeric(combat_BLBR_Beta_adjusted[y,BA7])))# across all BA7 samples
Covariance_ba10<-sapply(1:nrow(combat_BLBR_Beta_adjusted), function(y) Variation(as.numeric(combat_BLBR_Beta_adjusted[y,BA10])))# across all BA10 samples
Covariance_ba20<-sapply(1:nrow(combat_BLBR_Beta_adjusted), function(y) Variation(as.numeric(combat_BLBR_Beta_adjusted[y,BA20])))# across all BA20 samples

var_tissues<-data.frame(CpG=rownames(combat_BLBR_Beta_adjusted),
                        Var_blood=Covariance_blood,
                        Var_allbrain=Covariance_allbrain,
                        Var_BA7=Covariance_ba7,
                        Var_BA10=Covariance_ba10,
                        Var_BA20=Covariance_ba20)

save(var_tissues,file="CpG_Varbility_tissue.RData")



##### Sarah plot
load("CpG_Varbility_tissue.RData")

fit.data<-var_tissues


fit<-lm(Var_allbrain~0+Var_blood, data=fit.data)

fit.data<-fit.data[order(fit.data$Var_blood, decreasing=FALSE),]
head(fit.data)
fit.data$y1<-fit.data$Var_blood-(0.2)*fit.data$Var_blood
fit.data$y2<-fit.data$Var_blood+(0.2)*fit.data$Var_blood
head(fit.data)
fit.data$ind<-NA

fit.data

head(fit.data)
fit.data[, c(1, 6, 9)]
sum(is.na(fit.data$ind)) #36,113

plot(
  Var_allbrain ~ Var_blood,
  data = fit.data,
  col  = "black",    xlab="Brain Varibility", ylab="Blood Varibility", pch=20, cex=0.3)
lines(x=fit.data$Var_blood, y=fit.data$y1, col="grey")
lines(x=fit.data$Var_blood, y=fit.data$y2, col="grey")

ggplot()+geom_point(aes(-log(Var_blood), -log(Var_allbrain)), fit.data, shape=19, alpha=0.1)+theme_bw()+
  geom_line(aes(-log(Var_blood), -log(y1)), fit.data)+
  geom_line(aes(-log(Var_blood), -log(y2)), fit.data)
  
## varbility density plot
var_tissues_melt<-melt(var_tissues)

ggplot(var_tissues_melt, aes(log(value), color=variable))+geom_density()+theme_bw()
ggplot(var_tissues_melt, aes(variable,log(value), color=variable))+geom_boxplot()+theme_bw()

