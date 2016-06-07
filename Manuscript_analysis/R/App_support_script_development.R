load("cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     64
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta
load("SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta_adjusted[which(!(rownames(combat_BLBR_Beta_adjusted)%in%rownames(SnpatCpG))),]

correlations_BLBR<-read.csv("~/Documents/Blood_Brain/Python/correlation_jan22.csv", sep="\t")                       
colnames(correlations_BLBR)<-c("CpG","BRAIN7","BRAIN10","BRAIN20")
correlations_BLBR$CpG<-rownames(combat_BLBR_Beta)
correlations_BLBR<-correlations_BLBR[which(!(correlations_BLBR$CpG%in%rownames(SnpatCpG))),]

correlations_BLBR$MeanAllBrain<-rowMeans(correlations_BLBR[,2:4])
correlations_BLBR$sdMeanAllBrain<-sapply(1:nrow(correlations_BLBR), function(x) sd(correlations_BLBR[x,2:4]))

# percentile
pos<-quantile(correlations_BLBR$MeanAllBrain[which(correlations_BLBR$MeanAllBrain>=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
neg<-quantile(correlations_BLBR$MeanAllBrain[which(correlations_BLBR$MeanAllBrain<=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))

pos_percential<-unlist(sapply(1:nrow(correlations_BLBR), function(x) {corr<-correlations_BLBR$MeanAllBrain[x]
                                                                      if(corr<=0){"-"}else{
                                                                        if(corr<=pos[5]){"<50%"}else{
                                                                 if(corr<=pos[7]){"50-75%"}else{
                                                                   if(corr<=pos[9]){"75-90%"}else{"90%"}}}}
                                                               }))

neg_percential<-unlist(sapply(1:nrow(correlations_BLBR), function(x) {corr<-correlations_BLBR$MeanAllBrain[x]
                                                                      if(corr>=0){"-"}else{
                                                                      if(corr>=neg[5]){"<50%"}else{
                                                                        if(corr>=neg[2]){"50-75%"}else{
                                                                          if(corr>=neg[1]){"75-90%"}else{"90%"}}}}
                                                                      }))

correlations_BLBR$Pos_percentile<-pos_percential
correlations_BLBR$neg_percential<-neg_percential



meta<-read.csv("~/Documents/Blood_Brain/SampleInfo.csv")
fix_names<-unlist(lapply(meta$X, function(x) gsub(" ","", x , fixed=TRUE)))
meta$X<-fix_names
meta<-meta[which(meta$X%in%colnames(combat_BLBR_Beta)),]
meta<-meta[match(colnames(combat_BLBR_Beta), meta$X),]

load("Gene_CpG_Relations_updatejune2015.RData")
Gene_CpG_Relations_minimal<-Gene_CpG_Relations_update[,c(1:4,6)]
Gene_CpG_Relations_minimal<-Gene_CpG_Relations_minimal[!duplicated(Gene_CpG_Relations_minimal),]
CpG_gene<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs),]
## summarize the genes into a list
notdupli<-CpG_gene[which(!(CpG_gene$Probe_ID%in%CpG_gene$Probe_ID[duplicated(CpG_gene$Probe_ID)])),]
notdupli$gene<-as.character(notdupli$gene)
dupli<-CpG_gene[which(CpG_gene$Probe_ID%in%CpG_gene$Probe_ID[duplicated(CpG_gene$Probe_ID)]),]
if(nrow(dupli)>0){
  dupli<-split(dupli, dupli$Probe_ID)
  dupli<-lapply(1:length(dupli), function(x){y<-dupli[[x]][1,1:3]
                                             y$gene<-paste(as.character(dupli[[x]]$gene), sep=",",collapse=", ")
                                             y$region<-paste(as.character(dupli[[x]]$region), sep=",",collapse=", ")
                                             y})
  dupli<-do.call(rbind, dupli)
  CpG_gene<-rbind(notdupli, dupli)}else{}


# varibility
load("CpG_Varbility_tissue.RData")


save(combat_BLBR_Beta_adjusted,meta,correlations_BLBR,Gene_CpG_Relations_minimal,var_tissues,Gene_CpG_Relations_update, file="BLBR_app_Objects.RData")




##### gene selection
CpG_list<-function(gene_list, CpGs){
  gene_list<-as.character(gene_list)
  Genes<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$gene==gene_list | Gene_CpG_Relations_minimal$Probe_ID==CpGs),]
  
  Genes_onarray<-unique(Genes$gene)
  CpGs<-unique(c(as.character(Genes$Probe_ID), CpGs)) 
  withcorrelation<-nrow(correlations_BLBR[which(correlations_BLBR$CpG%in%CpGs),])
  paste(length(CpGs), "CpGs associated with ", length(Genes_onarray), "genes too look at (",
        withcorrelation, "CpGs have correlation data between blood and brain.", sep=" ")
}
CpG_list("AADACL4","cg00000109")


### make CpG list for other inputs
CpG_list_forplottable<-function(gene_list, CpGs){
  Genes<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$gene==gene_list | Gene_CpG_Relations_minimal$Probe_ID==CpGs),]
  
  Genes_onarray<-unique(Genes$gene)
  CpGs<-unique(c(as.character(Genes$Probe_ID), CpGs)) 
  CpGs
}
CpGs<-CpG_list_forplottable("AADACL4","cg00000109")


CpGs<-correlations_BLBR$CpG[1:input$CpGnum]








### summary table
summTable<-function(CpGs){
  ## associated Gene information

  ##correlation information
  correlations<-correlations_BLBR[which(correlations_BLBR$CpG%in%CpGs),]
  CpG_gene_correlations<-merge(CpG_gene, correlations, by.x="Probe_ID", by.y="CpG")
  CpG_gene_correlations$MeanAllBrain<-rowMeans(CpG_gene_correlations[,6:8])
  CpG_gene_correlations$sdMeanAllBrain<-sapply(1:nrow(CpG_gene_correlations), function(x) sd(CpG_gene_correlations[x,6:8]))
  
  ##varbility information
  varibility<-var_tissues[which(var_tissues$CpG%in%CpGs),]
  CpG_gene_correlations_varibility<-merge(CpG_gene_correlations, var_tissues, by.x="Probe_ID", by.y="CpG")
  CpG_gene_correlations_varibility<-CpG_gene_correlations_varibility[order(CpG_gene_correlations_varibility$Probe_ID),]
  CpG_gene_correlations_varibility[6:15]<-lapply(CpG_gene_correlations_varibility[6:15], function(x) round(x, digits=2))
  
  CpG_gene_correlations_varibility<-as.data.frame(lapply(CpG_gene_correlations_varibility, as.character))
  colnames(CpG_gene_correlations_varibility)<-c("CpG ID","Genomic Coordinate (build 37)","Chromosome (build 37)",
                                                "Associated Genes", "CpG in Feature of Gene, respectively","Correlation Blood-BA7",
                                                "Correlation Blood-BA10","Correlation Blood-BA20","Mean Correlation All Brain",
                                                "Standard Deviation Mean Correlation All Brain", "Varibility in Blood",
                                                "Varibility in All Brain","Varibility in BA7","Varibility in BA10","Varibility in BA20")
  CpG_gene_correlations_varibility
}


summTable(correlations_BLBR$CpG[1:10])



#################### plot comethylation

Comethylation_Plot<-function(correlations_BLBR, Betas, CpG_Hit_List, MaxCpGNum){
  hits_BLBR<-correlations_BLBR[which(correlations_BLBR$CpG%in%CpG_Hit_List),]
  if(nrow(hits_BLBR)>MaxCpGNum){hits_BLBR<-hits_BLBR[1:MaxCpGNum,]}else{}
  
  nrow(hits_BLBR)
  hits_BLBR_PBMC_correlation_melt<-melt(hits_BLBR)
  
  meta$SubNumber<-as.factor(meta$SubjectNumC)
  levels(meta$SubNumber)<-c(1:16)
  BLBR_Beta<-as.data.frame(Betas[which(rownames(Betas)%in%hits_BLBR$CpG),])
  BLBR_Beta$CpG<-rownames(BLBR_Beta)
  BLBR_Beta<-melt(BLBR_Beta)
  BLBR_Beta<-merge(BLBR_Beta, meta, by.x="variable", by.y="X")
  
  # Comparison highlight
  #hits_BLBR$max<-sapply(1:nrow(hits_BLBR), function(x) {
  #as.numeric(hits_BLBR[,c(2:4)][x,which(abs(hits_BLBR[x,c(2:4)])==max(abs(hits_BLBR[x,c(2:4)])))][1])})
  #hits_BLBR$max_abs<-abs(hits_BLBR$max)
  #BLBR_Beta$line<-sapply(1:nrow(BLBR_Beta), function(x) {
  #hit<-hits_BLBR[which(hits_BLBR$CpG==BLBR_Beta$CpG[x]),]
  #region_max_corr<-colnames(hit)[c(2:4)][which(hit[1,c(2:4)]==hit$max[1])]
  #if(as.character(BLBR_Beta$TissueType[x])==region_max_corr){"Max_Correlation"}else{"Less Correlated"}})
  #invisible(sapply(1:nrow(BLBR_Beta), function(x) if(BLBR_Beta$TissueType[x]=="PBMC"){
  #BLBR_Beta$line[x]<<-"Max_Correlation"}else{BLBR_Beta$line[x]<<-BLBR_Beta$line[x]}))
  
  ggplot()+geom_line(aes(SubNumber, value,group=TissueType,color=TissueType),#,alpha=line
                     BLBR_Beta, size=1.5)+
    theme_bw()+facet_wrap(~CpG)+
    #geom_text(data=hits_BLBR, aes(x=4, y=0.1, label=round(max, 2)),
    #colour="black", inherit.aes=FALSE, parse=FALSE)+
    #scale_alpha_manual(values = c(0.25, 1))+
    scale_color_manual(values=c("#fb6a4a","#ef3b2c","#cb181d","cornflowerblue"))+
    ylim(0,1)}

Comethylation_Plot(correlations_BLBR,combat_BLBR_Beta_adjusted, CpGs, 25)







##### correlation and varibilit in context

CpGs<-c("cg00308130","cg00201133")



#Place in header of server
correlations_BLBR_densityplt<-correlations_BLBR
correlations_BLBR_densityplt$mean<-rowMeans(correlations_BLBR_densityplt[,2:4])
  
  ##correlation information
cor_density_plot<-function(CpGs){
  correlations<-correlations_BLBR_densityplt[which(correlations_BLBR_densityplt$CpG%in%CpGs),]
  
  ggplot()+geom_density(aes(mean), correlations_BLBR_densityplt, fill="grey")+theme_bw()+
    geom_vline(aes(xintercept=mean), correlations)}



#Place in header of server
var_tissues_densityplt<-var_tissues[,1:3]
var_tissues_densityplt<-melt(var_tissues_densityplt, id="CpG")

##varbility information
var_density_plot<-function(CpGs){
  varibility<-var_tissues_densityplt[which(var_tissues_densityplt$CpG%in%CpGs),]
  
  ggplot()+geom_density(aes(value, fill=variable), var_tissues_densityplt)+theme_bw()+
    geom_vline(aes(xintercept=value), varibility)+facet_wrap(~variable, ncol=1)+
    theme(text = element_text(size=15))+xlab("CpG Varibility (Reference Range)")+
    ylab("Density")+scale_fill_manual(values=c("cornflowerblue","#fb6a4a"))}


  ##varbility information
  varibility<-var_tissues[which(var_tissues$CpG%in%CpGs),]
  CpG_gene_correlations_varibility<-merge(CpG_gene_correlations, var_tissues, by.x="Probe_ID", by.y="CpG")
  CpG_gene_correlations_varibility<-CpG_gene_correlations_varibility[order(CpG_gene_correlations_varibility$Probe_ID),]
  colnames(CpG_gene_correlations_varibility)<-c("CpG ID","Genomic Coordinate (build 37)","Chromosome (build 37)",
                                                "Associated Genes", "CpG in Feature of Gene, respectively","Correlation Blood-BA7",
                                                "Correlation Blood-BA10","Correlation Blood-BA20","Mean Correlation All Brain",
                                                "Standard Deviation Mean Correlation All Brain", "Varibility in Blood",
                                                "Varibility in All Brain","Varibility in BA7","Varibility in BA10","Varibility in BA20")
  CpG_gene_correlations_varibility
}