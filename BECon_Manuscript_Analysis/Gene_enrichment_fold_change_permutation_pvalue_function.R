




##### Permutation P value and fold change plot


CGI_Gene_permutation_enrichment<-function(CpG_list,background.probes, permutation_number, plotmin, plotmax){
  
  # GENE BODY
  Genes_correlated_CpGs<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%CpG_list),]
  Genes_correlated_CpGs<-Genes_correlated_CpGs[!duplicated(Genes_correlated_CpGs),]
  Genes_correlated_CpGs<-Genes_correlated_CpGs[!duplicated(Genes_correlated_CpGs[,c(1,4)]),]#remove duplicate CpG to gene associations
  
  Gene_hits_regionMeans<-tapply(Genes_correlated_CpGs$Probe_ID, Genes_correlated_CpGs$region, length)
  Gene_hits_regionMeans[is.na(Gene_hits_regionMeans)]<-0
  Gene_hits_regionMeans<-data.frame(Probe_Count=as.numeric(Gene_hits_regionMeans), Region=names(Gene_hits_regionMeans))
  
  ## Boot strapping (to see if hits more in feature than expected)
  bootstrap_genes<-lapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    Gene_rnd<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rnd_CpGs),]
    Gene_rnd<-Gene_rnd[!duplicated(Gene_rnd),]
    Gene_rnd<-Gene_rnd[!duplicated(Gene_rnd[,c(1,4)]),]#remove duplicate CpG to gene associations
    Gene_rnd_regionMeans<-tapply(Gene_rnd$Probe_ID, Gene_rnd$region, length)
    Gene_rnd_regionMeans[is.na(Gene_rnd_regionMeans)]<-0
    Gene_rnd_regionMeans
  })
  bootstrap_genes<-do.call(rbind, bootstrap_genes)
  
  print("FDR Adjusted Permutation P values for enrichment and depletion")
  Region_results<-sapply(1:nrow(Gene_hits_regionMeans), function(x){
    real_CpG_in_region<-Gene_hits_regionMeans$Probe_Count[x]
    Adjusted_enrich_p<-p.adjust(length(which(bootstrap_genes[,x]>=real_CpG_in_region))/permutation_number, method="fdr", n=nrow(Gene_hits_regionMeans))
    Adjusted_depletion_p<-p.adjust(length(which(bootstrap_genes[,x]<=real_CpG_in_region))/permutation_number, method="fdr", n=nrow(Gene_hits_regionMeans))
    print(paste("Enrichment: ", Adjusted_enrich_p, "; Depletion ", Adjusted_depletion_p, "; Feature: ", Gene_hits_regionMeans$Region[x],sep=""))
  })
  
  
  ## CGI
  CGI<-annotation[,c(58,49,50)]
  
  Resort_hits<-CGI[which(CGI$CpG%in%CpG_list),] #2138 CpGs, 1576 CGIs
  Resort_hits_featureMeans<-tapply(Resort_hits$CpG, Resort_hits$RELATION_TO_UCSC_CPG_ISLAND, length)
  Resort_hits_featureMeans[is.na(Resort_hits_featureMeans)]<-0
  Resort_hits_featureMeans<-data.frame(Probe_Count=as.numeric(Resort_hits_featureMeans),
                                       Feature=names(Resort_hits_featureMeans))
  levels(Resort_hits_featureMeans$Feature)<-c("None","Island","N_Shelf","N_Shore","S_Shelf","S_Shore")
  
  ## Boot strapping (to see if hits more in feature than expected)
  bootstrap_CGI<-lapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    Resort_rnd<-CGI[which(CGI$CpG%in%rnd_CpGs),]
    Resort_rnd_featureMeans<-tapply(Resort_rnd$CpG, Resort_rnd$RELATION_TO_UCSC_CPG_ISLAND, length)
    Resort_rnd_featureMeans[is.na(Resort_rnd_featureMeans)]<-0
    Resort_rnd_featureMeans
  })
  bootstrap_CGI<-do.call(rbind,bootstrap_CGI)
  
  
  CGI_results<-sapply(1:nrow(Resort_hits_featureMeans), function(x){
    real_CpG_in_region<-Resort_hits_featureMeans$Probe_Count[x]
    Adjusted_enrich_p<-p.adjust(length(which(bootstrap_CGI[,x]>=real_CpG_in_region))/permutation_number, method="fdr", n=nrow(Resort_hits_featureMeans))
    Adjusted_depletion_p<-p.adjust(length(which(bootstrap_CGI[,x]<=real_CpG_in_region))/permutation_number, method="fdr", n=nrow(Resort_hits_featureMeans))
    print(paste("Enrichment: ", Adjusted_enrich_p, "; Depletion: ", Adjusted_depletion_p, "; Feature: ", Resort_hits_featureMeans$Feature[x], sep=""))
  })
  
  c(Region_results,CGI_results)
  
  ### fold change plot
  
  ## fold changes
  colnames(Gene_hits_regionMeans)[2]<-"Feature"
  Gene_hits_regionMeans$Data<-"Gene Features"
  Resort_hits_featureMeans$Data<-"CpG Island Features"
  real<-rbind(Gene_hits_regionMeans,Resort_hits_featureMeans)
  print(real)
  
  boot<-cbind(bootstrap_genes,bootstrap_CGI)
  
  real$Fold_change<-sapply(1:nrow(real), function(x) mean(foldchange(real$Probe_Count[x], boot[,x])))
  real$Fold_change[is.infinite(real$Fold_change)]<-NA
  
  se <- function(x) sd(x)/sqrt(length(x))
  
  real$Fold_change_se<-sapply(1:nrow(real), function(x) se(foldchange(real$Probe_Count[x], boot[,x])))
  real$Fold_change_se[is.infinite(real$Fold_change_se)]<-NA
  
  real$Feature<-factor(real$Feature, levels=c("S_Shelf","S_Shore","Island","N_Shore","N_Shelf","promoter","intragenic","three_plus","intergenic","None"))
  
  levels(real$Feature)<-c("S. Shelf","S. Shore","Island","N. Shore","N. Shelf","Promoter","Intragenic","Three Prime","Intergenic","None")
  
  
  #set colors
  myColors <- c("#74add1","#feb24c","#66bd63","#feb24c","#74add1","#377eb8","#e41a1c","#984ea3","#999999","#999999")
  names(myColors) <- levels(real$Feature)
  fillscale <- scale_fill_manual(guide=F,values = myColors, drop = FALSE)
  
  
  
  print(ggplot(real, aes(Feature, Fold_change, fill=Feature))+
          geom_bar(position=position_dodge(width=0.9),stat="identity", color="black")+theme_bw()+
          fillscale+
          geom_errorbar(aes(ymax = Fold_change+Fold_change_se, ymin=Fold_change-Fold_change_se),width=0.25,position=position_dodge(width=0.9))+
          ylab("Fold Change")+xlab("")+ylim(plotmin,plotmax)+
          theme(legend.position="none",
                axis.text = element_text(size =12, color="black"),
                axis.title = element_text(size =15),
                legend.text = element_text(size =14),
                legend.title = element_text(size =12),
                strip.text.x = element_text(size = 12))+facet_wrap(~Data, scale="free_x"))
  
}
