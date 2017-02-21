library(GEOquery)
library(RCurl)
library(GEOmetadb)
library(reshape)
library(ggplot2)
setwd("/big_data/redgar/BECon/")


###############
## Get Brain Gene Expression Data
###############

getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
dbListFields(con, "gsm")

meta<-dbGetQuery(con, "select series_id,title,gpl,description,gsm,source_name_ch1,characteristics_ch1 from gsm")

# limit to two one expression
meta<-meta[which(meta$gpl%in%"GPL570"),] ## 314 healthy brain samples
brains<-meta[grep("Brain|brain|cortex|Cortex|striatum|hippocampus|cerebellum|frontal|Frontal|DLPFC",meta$source_name_ch1),]
brains$Tissue<-brains$source_name_ch1

# rename tissues to simplier name
brains$Tissue[grep("Brain|brain", brains$source_name_ch1)]<-"brain"
brains$Tissue[grep("Cortex|cortex", brains$source_name_ch1)]<-"cortex"
#BA10
brains$Tissue[grep("DLPFC", brains$source_name_ch1)]<-"prefrontal cortex"
#BA10
brains$Tissue[grep("frontal cortex|Frontal cortex|Frontal Cortex|frontal_cortex|Sporadic-frontal|Progranulin-frontal|Normal-frontal", brains$source_name_ch1)]<-"frontal cortex"
#BA20
brains$Tissue[grep("Temporal cortex|temporal cortex", brains$source_name_ch1)]<-"temporal cortex"
brains$Tissue[grep("entorhinal cortex|Entorhinal cortex", brains$source_name_ch1)]<-"entorhinal cortex"
brains$Tissue[grep("anterior cingulate cortex|Cingulate|cingulate", brains$source_name_ch1)]<-"cingulate cortex"
brains$Tissue[grep("amygdala|Amygdala", brains$source_name_ch1)]<-"amygdala"
brains$Tissue[grep("ganglia", brains$source_name_ch1)]<-"ganglia"
brains$Tissue[grep("striatum", brains$source_name_ch1)]<-"striatum"
brains$Tissue[grep("hippocampus|Hippocampus", brains$source_name_ch1)]<-"hippocampus"
#BA20
brains$Tissue[grep("temporal lobe of the cerebral cortex", brains$source_name_ch1)]<-"temporal lobe of the cerebral cortex"
brains$Tissue[grep("cerebellum|Cerebellum", brains$source_name_ch1)]<-"cerebellum"
#BA7
brains$Tissue[grep("parietal lobe of the cerebral cortex", brains$source_name_ch1)]<-"parietal lobe of the cerebral cortex"
#BA10
brains$Tissue[grep("frontal lobe of the cerebral cortex", brains$source_name_ch1)]<-"frontal lobe of the cerebral cortex"
brains$Tissue[grep("occipital lobe of the cerebral cortex", brains$source_name_ch1)]<-"occipital lobe of the cerebral cortex"
brains$Tissue[grep("cerebral cortex (region unknown)", brains$source_name_ch1)]<-"cerebral cortex (region unknown)"
brains$Tissue[grep("BA10", brains$source_name_ch1)]<-"BA10"
brains$Tissue[grep("gyrus|Gyrus", brains$source_name_ch1)]<-"gyrus"

# remove fetal cancers etc
brains<-brains[-grep("fetal|Fetal|Whole brain|gestational|rhesus|chimpanzee|glioma|brain,|cancer|meningioma|lesion|tumor|white_matter|Sporadic|carcinoma|pooled RNA|Reference RNA|astrocytoma|hemorrhage", brains$source_name_ch1),]
brains<-brains[-grep("fetal|Fetal|Whole brain|gestational|rhesus|chimpanzee|glioma|brain,|cancer|meningioma|lesion|tumor|white_matter|Sporadic|carcinoma|pooled RNA|Reference RNA|astrocytoma|hemorrhage", brains$description),]

unique(brains$Tissue)

#Since don't have much BA7 better strategy is to just pool all healthy brains to get a larger meta signal
brains$healthy<-"Unknown"
brains$healthy[grep("Control|control|healthy|Healthy", brains$description)]<-"Healthy"
brains$healthy[grep("Control|control|healthy|Healthy", brains$source_name_ch1)]<-"Healthy"
brain_healthy<-brains[which(brains$healthy=="Healthy"),]


### get actual expression data
GPL570_brains<-subset(brain_healthy, gpl=="GPL570")#174
GSE<-unique(GPL570_brains$series_id)
GEO<-lapply(GSE, function(x) as.data.frame(exprs(getGEO(x)[[1]])))
names(GEO)<-GSE
save(GEO, file="GEO_GPL570.RData") # 9 series; 54675 probes


###### Summarize the expression by each gene taking the mean of multiple probes
load("GEO_GPL570.RData")

GEO_GPL570<-GEO[[1]] # 1 series
# load the array annotation file (From GEO aswell)
GPL570_annotation<-read.csv("GPL570_GEO_annotation.csv", sep="\t")
GEO_GPL570$ID<-rownames(GEO_GPL570)
GEO_GPL570<-merge(GPL570_annotation[,c(1,11)], GEO_GPL570, by="ID")
#summarize to 21049 genes
GEO_GPL570_gene_sum<-lapply(3:ncol(GEO_GPL570), function(x) tapply(GEO_GPL570[,x], as.character(GEO_GPL570$Gene.Symbol), mean, na.rm=T))
GEO_GPL570_gene_sum<-as.data.frame(do.call(cbind, GEO_GPL570_gene_sum))
colnames(GEO_GPL570_gene_sum)<-colnames(GEO_GPL570)[3:ncol(GEO_GPL570)]
GEO_GPL570_gene_sum<-GEO_GPL570_gene_sum[-1,]#remove the no gene row

GEO_GPL570_gene_sum$gene<-rownames(GEO_GPL570_gene_sum)
Brain_series<-"GSE17612"
   
brain<-GEO_GPL570_gene_sum
brain_melt<-melt(brain, id="gene")
brain_melt$value<-as.numeric(as.character(brain_melt$value))
brain_melt$series<-"GSE17612"
ggplot(brain_melt, aes(log2(value), color=series))+geom_density()+theme_bw()+
  scale_color_manual(values=c("red","grey","grey","red",rep("grey",19)))

brain_melt<-brain_melt
brain_melt$tissue<-"brain"

brain_expression<-brain
save(brain_expression, file="Brain_gene_expression.RData")





###############
## Get Blood Gene Expression Data
###############
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
dbListFields(con, "gsm")

meta<-dbGetQuery(con, "select series_id,title,gpl,description,gsm,source_name_ch1,characteristics_ch1 from gsm")
# limit to two arrays
meta<-meta[which(meta$gpl%in%"GPL570"),]

blood<-meta[grep("WB|wb|blood|Blood",meta$source_name_ch1),]
blood$Tissue<-blood$source_name_ch1

#remove certain samples
blood<-blood[-grep("peripheral|Peripheral|leukocytes|Leukocytes|skin|MS|Bcells", blood$source_name_ch1),]
# rename to simplier name
blood$Tissue<-"blood"

# Healthy control filtering
blood$healthy<-"unknown"
blood$healthy[grep("Control|control|healthy|Healthy", blood$description)]<-"Healthy"
blood$healthy[grep("Control|control|healthy|Healthy", blood$source_name_ch1)]<-"Healthy"
blood$healthy[grep("Control|control|healthy|Healthy", blood$characteristics_ch1)]<-"Healthy"
blood_healthy<-blood[which(blood$healthy=="Healthy"),]


### get actual expression data
GPL570_blood<-subset(blood_healthy, gpl=="GPL570")#817
GPL6947_blood<-subset(blood_healthy, gpl=="GPL6947")#838


GSE<-unique(GPL570_blood$series_id)
GSE<-GSE[c(1:6,8:16,18:20,22,24:35)]
GSE<-GSE[c(1:5)]#mini

GEO<-lapply(GSE, function(x) as.data.frame(exprs(getGEO(x)[[1]])))
names(GEO)<-GSE
save(GEO, file="GEO_GPL570_blood.RData") # x series; x probes


### need more GPL570 blood so look for samples missed on first pass
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
dbListFields(con, "gsm")
x<-dbGetQuery(con, "select title,description,series_id,gsm,source_name_ch1,characteristics_ch1 from gsm where gpl='GPL570'")
meta<-x

blood<-meta[grep("WB|wb|blood|Blood",meta$source_name_ch1),]
blood$Tissue<-blood$source_name_ch1

#remove certain samples
blood<-blood[-grep("pheripheral|blood vessel|CD19+|CD1c+|Mixed sample|white blood|Tcells|White blood|MCF-7|Rhesus|CD8|lymphocytes|MCF-7|B cells|marrrow|Periferical|periperal|microvascular|Bone Marrow|newborn|Newborn|Natural Killer|Pluripotent|T cells|Cord|cord|umbilical|peripheral|Peripheral|leukocytes|Leukocytes|skin|MS|Bcells|monkey|macaque|rhesus|Diluted|monocyte", blood$source_name_ch1),]
# rename to simplier name
blood$Tissue<-"blood"

# Healthy control filtering
blood$healthy<-"unknown"
blood$healthy[grep("Control|control|healthy|Healthy", blood$description)]<-"Healthy"
blood$healthy[grep("Control|control|healthy|Healthy", blood$source_name_ch1)]<-"Healthy"
blood$healthy[grep("Control|control|healthy|Healthy", blood$characteristics_ch1)]<-"Healthy"
blood_healthy<-blood[which(blood$healthy=="Healthy"),]#687   8

save(blood_healthy, file="healthy_blood_expression.RData")

#GSE61635 and GSE55201

GSE<-c("GSE61635","GSE55201","GSE37171")
GEO<-lapply(GSE, function(x) as.data.frame(exprs(getGEO(x)[[1]])))
names(GEO)<-GSE
save(GEO, file="GEO_GPL570_bloodpart2.RData") # x series; x probes 


## BA 10 study is on GPL570 so just need blood from 570
load("healthy_blood_expression.RData")
load("GEO_GPL570_blood.RData")
GEO1<-GEO
load("GEO_GPL570_bloodpart2.RData")
GEO2<-GEO # one dataset is missing 62 genes

GEO1[6:8]<-GEO2
GEO<-lapply(1:length(GEO1), function(x) GEO1[[x]][which(rownames(GEO1[[x]])%in%rownames(GEO1[[8]])),])
GEO_GPL570<-GEO # 5 series for blood
GEO_GPL570<-do.call(cbind, GEO_GPL570)
GEO_GPL570<-GEO_GPL570[,which(colnames(GEO_GPL570)%in%blood_healthy$gsm)]# 54613   316 samples
blood_healthy<-blood_healthy[which(blood_healthy$gsm%in%colnames(GEO_GPL570)),]



###### Summarize the expression by each gene taking the mean of multiple probes

## add annotation to data
GPL570_annotation<-read.csv("GPL570_GEO_annotation.csv", sep="\t")
GEO_GPL570$ID<-rownames(GEO_GPL570)
GEO_GPL570<-merge(GPL570_annotation[,c(1,11)], GEO_GPL570, by="ID")
#summarize to 21049 genes
GEO_GPL570_gene_sum<-lapply(3:ncol(GEO_GPL570), function(x) tapply(GEO_GPL570[,x], as.character(GEO_GPL570$Gene.Symbol), mean, na.rm=T))
GEO_GPL570_gene_sum<-as.data.frame(do.call(cbind, GEO_GPL570_gene_sum))
colnames(GEO_GPL570_gene_sum)<-colnames(GEO_GPL570)[3:ncol(GEO_GPL570)]
GEO_GPL570_gene_sum<-GEO_GPL570_gene_sum[-1,]#remove the no gene row

GEO_GPL570_gene_sum$gene<-rownames(GEO_GPL570_gene_sum)
GEO_blood_genes<-GEO_GPL570_gene_sum
 
blood_healthy_plot<-blood_healthy[,c(3,4,7,8)]
   
blood<-GEO_GPL570_gene_sum
melt_blood<-melt(blood, id="gene")
melt_blood<-merge(melt_blood, blood_healthy_plot, by.x="variable", by.y="gsm")

ggplot(melt_blood, aes(log2(value), color=series_id))+geom_density()+theme_bw()+
  scale_color_manual(values=c("red","red","red","red","red","red","blue","blue","red","blue","blue",rep("grey",19)))


################ compare blood and brain
melt_blood$tissue<-"blood"
melt_blood<-melt_blood[,c(2,1,3,4,7)]
colnames(melt_blood)[4]<-"series"

mini_melt_both<-rbind(melt_blood, mini_melt_brain)

ggplot(mini_melt_both, aes(log2(value), color=series))+geom_density()+theme_bw()+scale_color_manual(values=c("black", "cornflowerblue", "green","grey","pink","blue","red"))
ggplot(mini_melt_both[which(mini_melt_both$value<15),], aes((value), color=series))+geom_density()+theme_bw()+scale_color_manual(values=c("black", "cornflowerblue", "green","grey","pink","blue","red"))

# some series had extermly different mean expression levels across all genes and very differetn distributions compared to the brain dataset available
# so it is likely in some data sets the data was transfromed in some way and therefore we selceted to blood studies with distributions most similar to the brain study
GEO_blood_genes<-GEO_blood_genes[, 1:316]
good_blood_GEO<-GEO_blood_genes[,which(blood_healthy$series_id%in%c("GSE37171","GSE61635"))]#21049    70

save(good_blood_GEO, file="good_blood_GEO_expressionGPL570.RData")







###############
## Informative Gene Expression in the collected blood and brain data
###############
BLBR_genes<-read.csv("Informative_CpGs_MultipleHitGenes_pos_neg_allchr.csv")# 10 CpGs per gene (239 genes)


#### Brain
load("Brain_gene_expression.RData")
load("Gene_annotation/Gene_CpG_Relations_updatejune2015.RData")

BLBR_genes<-BLBR_genes$x

# make gene name the first column
brain_expression<-brain_expression[,c(52, 1:(ncol(brain_expression)-1))]

#189/237 genes have expression data
GEO_brain_genes_Meth<-brain_expression[which(brain_expression$gene%in%BLBR_genes),]
# or 80% of genes
(nrow(GEO_brain_genes_Meth)/length(BLBR_genes))*100

## take a mean expression across all samples
BLBR_expression<-data.frame(mn=rowMeans(GEO_brain_genes_Meth[,2:ncol(GEO_brain_genes_Meth)]),
                            sd=sapply(1:nrow(GEO_brain_genes_Meth), function(x) sd(GEO_brain_genes_Meth[x,2:ncol(GEO_brain_genes_Meth)])),
                            Data="Blood Brain Correlated Genes")

## Build a random expected level of expression 
RND_expression<-lapply(1:100, function(x) {
  set.seed(x)
  Random_genes<-Gene_CpG_Relations_update$gene[sample(1:nrow(Gene_CpG_Relations_update),length(BLBR_genes))]
  GEO_rnd_genes_Meth<-brain_expression[which(brain_expression$gene%in%Random_genes),]
  print((nrow(GEO_rnd_genes_Meth)/length(Random_genes))*100)
  RND_expression<-data.frame(mn=rowMeans(GEO_rnd_genes_Meth[,2:ncol(GEO_rnd_genes_Meth)]),
                             sd=sapply(1:nrow(GEO_rnd_genes_Meth), function(x) sd(GEO_rnd_genes_Meth[x,2:ncol(GEO_rnd_genes_Meth)])),
                             Data="Random Genes")
  RND_expression})

RND_expression<-do.call(rbind,RND_expression)

plot<-data.frame(Mean=c(mean(RND_expression[,1]),mean(BLBR_expression[,1])),
                 SD=c(sd(RND_expression[,2]),sd(BLBR_expression[,2])),
                 Data=c("RND","BLBR"))

ggplot(plot, aes(Data, Mean, fill=Data))+geom_bar(position=position_dodge(width=0.9),stat="identity", color="black")+theme_bw()+
  geom_errorbar(aes(ymax = Mean + SD, ymin=Mean - SD),width=0.25,position=position_dodge(width=0.9))

Expression_compare_brain<-rbind(BLBR_expression, RND_expression)



### Blood
load("good_blood_GEO_expressionGPL570.RData")
load("Gene_annotation/Gene_CpG_Relations_updatejune2015.RData")

BLBR_genes<-BLBR_genes$x

# make gene name the first column
good_blood_GEO$gene<-rownames(good_blood_GEO)
good_blood_GEO<-good_blood_GEO[,c(71, 1:(ncol(good_blood_GEO)-1))]

#189/237 genes have expression data
GEO_blood_genes_Meth<-good_blood_GEO[which(good_blood_GEO$gene%in%BLBR_genes),]
# or 80% of genes
(nrow(GEO_blood_genes_Meth)/length(BLBR_genes))*100

BLBR_expression<-data.frame(mn=rowMeans(GEO_blood_genes_Meth[,2:ncol(GEO_blood_genes_Meth)]),
                            sd=sapply(1:nrow(GEO_blood_genes_Meth), function(x) sd(GEO_blood_genes_Meth[x,2:ncol(GEO_blood_genes_Meth)])),
                            Data="Blood Brain Correlated Genes")

## Build a random expected level of expression
RND_expression<-lapply(1:100, function(x) {
  set.seed(x)
  Random_genes<-Gene_CpG_Relations_update$gene[sample(1:nrow(Gene_CpG_Relations_update),length(BLBR_genes))]
  GEO_rnd_genes_Meth<-good_blood_GEO[which(good_blood_GEO$gene%in%Random_genes),]
  print((nrow(GEO_rnd_genes_Meth)/length(Random_genes))*100)
  RND_expression<-data.frame(mn=rowMeans(GEO_rnd_genes_Meth[,2:ncol(GEO_rnd_genes_Meth)]),
                             sd=sapply(1:nrow(GEO_rnd_genes_Meth), function(x) sd(GEO_rnd_genes_Meth[x,2:ncol(GEO_rnd_genes_Meth)])),
                             Data="Random Genes")
  RND_expression})

RND_expression<-do.call(rbind,RND_expression)

plot<-data.frame(Mean=c(mean(RND_expression[,1]),mean(BLBR_expression[,1])),
                 SD=c(sd(RND_expression[,2]),sd(BLBR_expression[,2])),
                 Data=c("RND","BLBR"))

ggplot(plot, aes(Data, Mean, fill=Data))+geom_bar(position=position_dodge(width=0.9),stat="identity", color="black")+theme_bw()+
  geom_errorbar(aes(ymax = Mean + SD, ymin=Mean - SD),width=0.25,position=position_dodge(width=0.9))


Expression_compare_blood<-rbind(BLBR_expression, RND_expression)



### Gene expression blood brain comparison
Expression_compare_blood$Tissue<-"Blood"
Expression_compare_brain$Tissue<-"Brain"
expression_compare<-rbind(Expression_compare_blood, Expression_compare_brain)

levels(expression_compare$Data)<-c("Blood Brain\nCorrelated Genes","Random Genes")

ggplot(expression_compare, aes(Data, log2(mn)))+geom_violin(fill="lightgrey", color="white")+
  geom_boxplot(aes(fill=Data),color="black", width=0.5, outlier.shape=19, outlier.size=1)+theme_bw()+facet_wrap(~Tissue)+
  ylab("Mean log2 Expression \n (Each Concordant Gene Across Gene Expression Samples)")+
  scale_fill_manual(values=c("cornflowerblue","#deebf7"))

Exp_blood_actual<-Expression_compare_blood$mn[which(Expression_compare_blood$Data=="Blood Brain Correlated Genes")]
Exp_blood_rnd<-Expression_compare_blood$mn[which(Expression_compare_blood$Data=="Random Genes")]

Exp_brain_actual<-Expression_compare_brain$mn[which(Expression_compare_brain$Data=="Blood Brain Correlated Genes")]
Exp_brain_rnd<-Expression_compare_brain$mn[which(Expression_compare_brain$Data=="Random Genes")]

wilcox.test(Exp_blood_actual, Exp_blood_rnd)
wilcox.test(Exp_brain_actual, Exp_brain_rnd)



