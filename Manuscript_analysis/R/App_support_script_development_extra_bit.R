# library(reshape)
# library(ggplot2)
# library(RColorBrewer)
# 
# setwd("/big_data/redgar/echelon_backup_2/Blood_Brain")
# load("BLBR_app_Objects.RData")
# 
# correlations_BLBR$Data<-"correlations"
# var_tissues$Data<-"varibility"
# 
# 
# colnames(correlations_BLBR)<-c("CpG","Correlation_BA7", "Correlation_BA10","Correlation_BA20",
#                               "MeanAllBrain","sdMeanAllBrain","Pos_percentile","neg_percential","Data")
# 
# cor_for_plot<-correlations_BLBR[,c(1:4,9)]
# 
# 
# ### Correlation percentile colors
# pos7<-quantile(cor_for_plot$Correlation_BA7[which(cor_for_plot$Correlation_BA7>=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
# pos10<-quantile(cor_for_plot$Correlation_BA10[which(cor_for_plot$Correlation_BA10>=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
# pos20<-quantile(cor_for_plot$Correlation_BA20[which(cor_for_plot$Correlation_BA20>=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
# 
# neg7<-quantile(cor_for_plot$Correlation_BA7[which(cor_for_plot$Correlation_BA7<=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
# neg10<-quantile(cor_for_plot$Correlation_BA10[which(cor_for_plot$Correlation_BA10<=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
# neg20<-quantile(cor_for_plot$Correlation_BA20[which(cor_for_plot$Correlation_BA20<=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
# 
#             
# 
# percential<-function(correlations, levelpos, levelneg) {sapply(correlations, function(corr){
#   pos<-levelpos
#   neg<-levelneg
#   if(corr>=0){
#     if(corr<=pos[5]){"<50% (positive)"}else{
#       if(corr<=pos[7]){"50-75% (positive)"}else{
#         if(corr<=pos[9]){"75-90% (positive)"}else{"90% (positive)"}}}}else{
#     if(corr>=neg[5]){"<50% (negative)"}else{
#       if(corr>=neg[2]){"50-75% (negative)"}else{
#         if(corr>=neg[1]){"75-90% (negative)"}else{"90% (negative)"}}}}})}
# 
# 
# cor_for_plot$percentile7<-percential(cor_for_plot$Correlation_BA7, pos7, neg7)
# cor_for_plot$percentile10<-percential(cor_for_plot$Correlation_BA10, pos10,neg10)
# cor_for_plot$percentile20<-percential(cor_for_plot$Correlation_BA20, pos20,neg20)
# 
# cor_for_plot_percentile<-melt(cor_for_plot[,c(1,6:8)], id="CpG")
# cor_for_plot_correlations<-melt(cor_for_plot[,c(1:5)], id=c("CpG","Data"))
# 
# cor_for_plot<-cbind(cor_for_plot_correlations, cor_for_plot_percentile)
# colnames(cor_for_plot)[7]<-"color"
# 
# cor_for_plot<-cor_for_plot[,c(1:4,7)]
# 
# 
# ### varibility threshold colors
# colnames(var_tissues)<-c("CpG","Varibility_Blood", "Var_allbrain","Varibility_BA7","Varibility_BA10","Varibility_BA20","Data")
# var_for_plot<-melt(var_tissues[,c(1,2,4:6,7)], id=c("CpG","Data"))
# var_for_plot$color<-sapply(1:nrow(var_for_plot), function(x) if(var_for_plot$value[x]>0.05){"variable"}else{"not variable"})
# 
# 
# #bind the data
# cor_var_for_plot<-rbind(cor_for_plot, var_for_plot)
# 
# cor_var_for_plot$color<-factor(cor_var_for_plot$color, levels=c("90% (positive)", "75-90% (positive)", "50-75% (positive)","<50% (positive)",
#                                                                 "<50% (negative)","50-75% (negative)","75-90% (negative)", "90% (negative)",
#                                                                 "variable","not variable"))
# 
# cor_var_for_plot$Data<-factor(cor_var_for_plot$Data,levels=c("varibility","correlations"))
# 
# #round for plotting
# cor_var_for_plot[,4]<-sapply(cor_var_for_plot[4], function(x) round(x, digits=2))
# 
# 
# 
# 
# 
# 
# ### add genomic information
# load("Gene_CpG_Relations_updatejune2015.RData")
# 
# Gene_CpG_Relations_minimal_BLBR<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%correlations_BLBR$CpG),]
# CpG_gene_melt<-melt(Gene_CpG_Relations_minimal_BLBR, id="Probe_ID")
# 
# CpG_gene_melt$variable<-factor(CpG_gene_melt$variable, levels=c("Chromosome_37","Coordinate_37","gene", "region"))
# levels(CpG_gene_melt$variable)<-c("Chromosome","Coordinate","Gene(s)", "Gene Region(s)")
# 
# colnames(CpG_gene_melt)<-c("CpG","Data","value")
# CpG_gene_melt$variable<-"NA"
# CpG_gene_melt$color<-"genomic"
# 
# CpG_gene_melt<-CpG_gene_melt[,c(1,2,4,3,5)]
# 
# 
# 
# CpG_gene_melt_combined<-rbind(CpG_gene_melt, cor_var_for_plot)
# CpG_gene_melt_combined$Data<-as.factor(CpG_gene_melt_combined$Data)
# 
# ## add cell composition
# load("Brain_Blood_changes_for_app.RData")
# Brain_Blood_changes_melt<-Brain_Blood_changes_melt[which(Brain_Blood_changes_melt$CpG%in%correlations_BLBR$CpG),]
# 
# CpG_gene_melt_combined<-rbind(CpG_gene_melt_combined, Brain_Blood_changes_melt)
# 
# 
# CpG_gene_melt_combined$Data<-factor(CpG_gene_melt_combined$Data,levels=c("Chromosome","Coordinate","Gene(s)","Gene Region(s)","varibility","correlations","Cell_composition"))
# CpG_gene_melt_combined$variable<-as.factor(CpG_gene_melt_combined$variable)
# 
# 
# 
# load("BLBR_app_Objects.RData")
# 
# # put things in a prettier order
# CpG_gene_melt_combined$color<-as.factor(CpG_gene_melt_combined$color)
# 
# CpG_gene_melt_combined$color<-factor(CpG_gene_melt_combined$color, 
#                                      levels=levels(CpG_gene_melt_combined$color)[c(9,10,11,1,2,3,4,5,6,7,8,18,17,16,19,14,13,12,15)])
# 
# 
# levels(CpG_gene_melt_combined$color)<-c("Genominc Info", "Variable","Not Variable",
#                                         "90% (Positive)","75-90% (Positive)",
#                                         "50-75% (Positive)","<50% (Positive)",
#                                         "<50% (Negative)","50-75% (Negative)",
#                                         "75-90% (Negative)","90% (Negative)",
#                                         "90% (Brain Cell Comp.)", 
#                                         "75-90% (Brain Cell Comp.)",
#                                         "50-75% (Brain Cell Comp.)",
#                                         "<50% (Brain Cell Comp.)", 
#                                         "90% (Blood Cell Comp.)",
#                                         "75-90% (Blood Cell Comp.)",
#                                         "50-75% (Blood Cell Comp.)",
#                                         "<50% (Blood Cell Comp.)")
# 
# save(combat_BLBR_Beta_adjusted,meta,correlations_BLBR,Gene_CpG_Relations_minimal,var_tissues,Gene_CpG_Relations_update,CpG_gene_melt_combined, file="BLBR_app_Objects.RData")
# 
# 
# 
# ## This should be in the server.R script
# 
# 
# # COLORS! define here so they are consistent between plots
# myColors <- c("white","#74c476","#c7e9c0",
#               "#fe9929","#fec44f",
#               "#fee391","#ffffe5",
#               "#f0f0f0","#bdbdbd",
#               "#969696","#737373",
#               "#4292c6","#6baed6",
#               "#9ecae1","#deebf7",
#               "#ef3b2c","#fb6a4a",
#               "#fc9272","#fee0d2")
# 
# names(myColors) <- levels(CpG_gene_melt_combined$color)
# fillscale <- scale_fill_manual(name = "Correlation or Cell Composition Percentile\n or Varibility Status",values = myColors, drop = FALSE)
# 
# 
# ##### PLot function
# Colorful_table_plot<-function(CpGs){
#   CpG_gene_melt_combined_mini<-CpG_gene_melt_combined[which(CpG_gene_melt_combined$CpG%in%CpGs),]
#   
#   levels(CpG_gene_melt_combined_mini$variable)<-c("Blood","Brain","BA10","BA20","BA7","","BA10","BA20","BA7","Blood")
#   levels(CpG_gene_melt_combined_mini$Data)<-c("Chromosome","Coordinate","Gene(s)","Gene Region(s)","Varibility","Correlation", "Cell Composition")
# 
#   CpG_gene_melt_combined_mini$CpG<-as.factor(CpG_gene_melt_combined_mini$CpG)
#   #CpG order
#   CpG_order<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs),]
#   CpG_order<-CpG_order[with(CpG_order, order(gene, Coordinate_37)), ]
#   
#   CpG_gene_melt_combined_mini$CpG<-factor(CpG_gene_melt_combined_mini$CpG, 
#                                           levels=CpG_order$Probe_ID)
#   
#   ggplot(CpG_gene_melt_combined_mini, aes(variable, CpG, fill = color)) +
#     geom_tile(color = "black",size=0.5) +
#     theme_gray(8)+facet_grid(.~Data, scales="free_x",space="free_x")+
#     geom_text(aes(label=value),color="black")+ #aes(label = round(value,2)), 
#     fillscale+xlab("")+ylab("")+
#     theme_minimal()+ guides(fill=guide_legend(ncol=2))+
#     theme(axis.text.x = element_text(size =12, color="black",margin=margin(-15,0,0,0)),
#           axis.ticks = element_blank(),
#           axis.title = element_text(size =15),
#           plot.margin=unit(c(15,0,10,0),"mm"),
#           legend.text = element_text(size =10),
#           legend.title = element_text(size =12),
#           strip.text.x = element_text(size = 12, face="bold"))
# }
# 
# 
# #test plot
# CpGs<-sample(CpG_gene_melt_combined$CpG, 5)
# Colorful_table_plot(CpGs)
# 
# 







############# Add all the info to the correlations_BLBR object then subset and melt in the server function
library(reshape)
library(ggplot2)
library(RColorBrewer)

setwd("/big_data/redgar/echelon_backup_2/Blood_Brain")
load("BLBR_app_Objects.RData")

correlations_BLBR$Data<-"correlations"
var_tissues$Data<-"varibility"


# colnames(correlations_BLBR)<-c("CpG","Correlation_BA7", "Correlation_BA10","Correlation_BA20",
#                                "MeanAllBrain","sdMeanAllBrain","Pos_percentile","neg_percential","Data")



### Correlation percentile colors
pos7<-quantile(correlations_BLBR$Correlation_BA7[which(correlations_BLBR$Correlation_BA7>=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
pos10<-quantile(correlations_BLBR$Correlation_BA10[which(correlations_BLBR$Correlation_BA10>=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
pos20<-quantile(correlations_BLBR$Correlation_BA20[which(correlations_BLBR$Correlation_BA20>=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))

neg7<-quantile(correlations_BLBR$Correlation_BA7[which(correlations_BLBR$Correlation_BA7<=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
neg10<-quantile(correlations_BLBR$Correlation_BA10[which(correlations_BLBR$Correlation_BA10<=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))
neg20<-quantile(correlations_BLBR$Correlation_BA20[which(correlations_BLBR$Correlation_BA20<=0)], c(0.1, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.8, 0.9))



percential<-function(correlations, levelpos, levelneg) {sapply(correlations, function(corr){
  pos<-levelpos
  neg<-levelneg
  if(corr>=0){
    if(corr<=pos[5]){"<50% (Positive)"}else{
      if(corr<=pos[7]){"50-75% (Positive)"}else{
        if(corr<=pos[9]){"75-90% (Positive)"}else{"90% (Positive)"}}}}else{
          if(corr>=neg[5]){"<50% (Negative)"}else{
            if(corr>=neg[2]){"50-75% (Negative)"}else{
              if(corr>=neg[1]){"75-90% (Negative)"}else{"90% (Negative)"}}}}})}


correlations_BLBR$percentile7<-percential(correlations_BLBR$Correlation_BA7, pos7, neg7)
correlations_BLBR$percentile10<-percential(correlations_BLBR$Correlation_BA10, pos10,neg10)
correlations_BLBR$percentile20<-percential(correlations_BLBR$Correlation_BA20, pos20,neg20)


### Varibility
#colnames(var_tissues)<-c("CpG","Varibility_Blood", "Var_allbrain","Varibility_BA7","Varibility_BA10","Varibility_BA20","Data")
# var_for_plot<-melt(var_tissues[,c(1,2,4:6,7)], id=c("CpG","Data"))

#variable yes/no function
var_y_n<-function(varibility) {if(varibility>0.05){"variable"}else{"not variable"}}

var_tissues$Var_blood_color<-sapply(var_tissues$Varibility_Blood, var_y_n)
var_tissues$Var_BA7_color<-sapply(var_tissues$Varibility_BA7, var_y_n)
var_tissues$Var_BA10_color<-sapply(var_tissues$Varibility_BA10, var_y_n)
var_tissues$Var_BA20_color<-sapply(var_tissues$Varibility_BA20, var_y_n)

## Cell composition

load("Brain_Blood_changes_for_app.RData")
# Brain_Blood_changes_melt<-Brain_Blood_changes_melt[which(Brain_Blood_changes_melt$CpG%in%correlations_BLBR$CpG),]

Brain_Blood_changes_variable<-cast(Brain_Blood_changes_melt, CpG~variable)

Brain_Blood_changes_color<-cast(Brain_Blood_changes_melt[,c(1,3,5)], CpG~variable)
colnames(Brain_Blood_changes_color)<-c("CpG","Blood_CellComp_Percentile","Brain_CellComp_Percentile")


Brain_Blood_changes<-merge(Brain_Blood_changes_variable, Brain_Blood_changes_color, by="CpG")

correlations_BLBR<-merge(correlations_BLBR, Brain_Blood_changes, by="CpG")




## rename levels
correlations_BLBR$percentile7<-as.factor(correlations_BLBR$percentile7)
correlations_BLBR$percentile7<-factor(correlations_BLBR$percentile7, 
                                      levels=c("90% (Positive)", "75-90% (Positive)", "50-75% (Positive)", "<50% (Positive)",
                                               "<50% (Negative)", "50-75% (Negative)", "75-90% (Negative)", "90% (Negative)"))
correlations_BLBR$percentile10<-as.factor(correlations_BLBR$percentile10)
correlations_BLBR$percentile10<-factor(correlations_BLBR$percentile10, 
                                      levels=c("90% (Positive)", "75-90% (Positive)", "50-75% (Positive)", "<50% (Positive)",
                                               "<50% (Negative)", "50-75% (Negative)", "75-90% (Negative)", "90% (Negative)"))
correlations_BLBR$percentile20<-as.factor(correlations_BLBR$percentile20)
correlations_BLBR$percentile20<-factor(correlations_BLBR$percentile20, 
                                      levels=c("90% (Positive)", "75-90% (Positive)", "50-75% (Positive)", "<50% (Positive)",
                                               "<50% (Negative)", "50-75% (Negative)", "75-90% (Negative)", "90% (Negative)"))



correlations_BLBR$Blood_CellComp_Percentile<-as.factor(correlations_BLBR$Blood_CellComp_Percentile)
levels(correlations_BLBR$Blood_CellComp_Percentile)<-c("50-75% (Blood Cell Comp.)","75-90% (Blood Cell Comp.)","90% (Blood Cell Comp.)","<50% (Blood Cell Comp.)")
correlations_BLBR$Blood_CellComp_Percentile<-factor(correlations_BLBR$Blood_CellComp_Percentile, levels=c("90% (Blood Cell Comp.)","75-90% (Blood Cell Comp.)",
                                                                                                          "50-75% (Blood Cell Comp.)","<50% (Blood Cell Comp.)"))

correlations_BLBR$Brain_CellComp_Percentile<-as.factor(correlations_BLBR$Brain_CellComp_Percentile)
levels(correlations_BLBR$Brain_CellComp_Percentile)<-c("50-75% (Brain Cell Comp.)","75-90% (Brain Cell Comp.)","90% (Brain Cell Comp.)","<50% (Brain Cell Comp.)")
correlations_BLBR$Brain_CellComp_Percentile<-factor(correlations_BLBR$Brain_CellComp_Percentile, levels=c("90% (Brain Cell Comp.)","75-90% (Brain Cell Comp.)",
                                                                                                          "50-75% (Brain Cell Comp.)","<50% (Brain Cell Comp.)"))

var_tissues$Var_BA7_color<-as.factor(var_tissues$Var_BA7_color)
levels(var_tissues$Var_BA7_color)<-c("Not Variable","Variable")
var_tissues$Var_BA10_color<-as.factor(var_tissues$Var_BA10_color)
levels(var_tissues$Var_BA10_color)<-c("Not Variable","Variable")
var_tissues$Var_BA20_color<-as.factor(var_tissues$Var_BA20_color)
levels(var_tissues$Var_BA20_color)<-c("Not Variable","Variable")
var_tissues$Var_blood_color<-as.factor(var_tissues$Var_blood_color)
levels(var_tissues$Var_blood_color)<-c("Not Variable","Variable")


########
# SAVE HERE
#######
save(combat_BLBR_Beta_adjusted,meta,correlations_BLBR,Gene_CpG_Relations_minimal,var_tissues, file="BLBR_app_Objects_beta.RData")




# COLORS! define here so they are consistent between plots
myColors <- c("white","#c7e9c0","#74c476",
              "#fe9929","#fec44f",
              "#fee391","#ffffe5",
              "#f0f0f0","#bdbdbd",
              "#969696","#737373",
              "#4292c6","#6baed6",
              "#9ecae1","#deebf7",
              "#ef3b2c","#fb6a4a",
              "#fc9272","#fee0d2")

color_possibilities<-c("Genomic Info", levels(var_tissues$Var_blood_color), levels(correlations_BLBR$percentile7), 
                       levels(correlations_BLBR$Blood_CellComp_Percentile), levels(correlations_BLBR$Brain_CellComp_Percentile))

names(myColors) <- color_possibilities
fillscale <- scale_fill_manual(name = "Correlation or Cell Composition Percentile\n or Varibility Status",values = myColors, drop = FALSE)





###### FUnction to plot
CpGs<-c("cg00000165","cg00000289","cg00000236")

Colorful_table_plot(CpGs)

Colorful_table_plot<-function(CpGs){
  color_table_plot_cor<-correlations_BLBR[which(correlations_BLBR$CpG%in%CpGs),c(1:4,10:16)]
  color_table_plot_var<-var_tissues[which(var_tissues$CpG%in%CpGs),c(1,2,4:6,8:11)]
  
  color_table_plot_cor_values<-melt(color_table_plot_cor[,c(1:4,8,9)])
  color_table_plot_var_values<-melt(color_table_plot_var[,c(1:5)])
  
  color_table_plot_cor_color<-melt(color_table_plot_cor[,c(1,5:7,10:11)], id="CpG")
  color_table_plot_var_color<-melt(color_table_plot_var[,c(1,6:9)], id="CpG")
  
  color_table_plot_values<-rbind(color_table_plot_cor_values, color_table_plot_var_values)
  color_table_plot_values$Data<-as.factor(color_table_plot_values$variable)
  levels(color_table_plot_values$Data)<-c("correlations","correlations","correlations","Cell_composition","Cell_composition","varibility","varibility","varibility","varibility")
  color_table_plot_color<-rbind(color_table_plot_cor_color, color_table_plot_var_color)
  color_table_plot_values$color<-color_table_plot_color$value
  color_table_plot_values<-color_table_plot_values[,c(1,4,2,3,5)]
  
  #round for plotting once leted in subset
  color_table_plot_values$value<-sapply(color_table_plot_values$value, function(x) round(x, digits=2)) 
  
  ## add genomic info
  ### rbind genomic information
  Gene_CpG_Relations_minimal_BLBR<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%color_table_plot_cor$CpG),]
  CpG_gene_melt<-melt(Gene_CpG_Relations_minimal_BLBR, id="Probe_ID")
  CpG_gene_melt$variable<-factor(CpG_gene_melt$variable, levels=c("Chromosome_37","Coordinate_37","gene", "region"))
  levels(CpG_gene_melt$variable)<-c("Chromosome","Coordinate","Gene(s)", "Gene Region(s)")
  colnames(CpG_gene_melt)<-c("CpG","Data","value")
  CpG_gene_melt$variable<-"NA"
  CpG_gene_melt$color<-"Genomic Info"
  CpG_gene_melt<-CpG_gene_melt[,c(1,2,4,3,5)]
  
  ## final data
  cor_var_for_plot<-rbind(color_table_plot_values, CpG_gene_melt)
  
  

  
  levels(cor_var_for_plot$variable)<-c("BA7","BA10","BA20","Blood","Brain","Blood","BA7","BA10","BA20","")
  levels(cor_var_for_plot$Data)<-c("Correlation","Cell Composition","Varibility","Chr","Coor","Gene(s)","Gene \nRegion(s)")
  
  cor_var_for_plot$Data<-factor(cor_var_for_plot$Data, levels=c("Chr","Coor","Gene(s)","Gene \nRegion(s)",
                                                                "Varibility","Correlation","Cell Composition"))
  
  cor_var_for_plot$CpG<-as.factor(cor_var_for_plot$CpG)
  #CpG order
  CpG_order<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs),]
  CpG_order<-CpG_order[with(CpG_order, order(gene, Coordinate_37)), ]
  
  cor_var_for_plot$CpG<-factor(cor_var_for_plot$CpG, 
                                          levels=CpG_order$Probe_ID)
  
  ggplot(cor_var_for_plot, aes(variable, CpG, fill = color)) +
    geom_tile(color = "black",size=0.5) +
    theme_gray(8)+facet_grid(.~Data, scales="free_x",space="free_x")+
    geom_text(aes(label=value),color="black")+ #aes(label = round(value,2)), 
    fillscale+xlab("")+ylab("")+
    theme_minimal()+ guides(fill=guide_legend(ncol=2))+
    theme(axis.text.x = element_text(size =12, color="black",margin=margin(-15,0,0,0), face="bold"),
          axis.text.y = element_text(size =12, color="black"),
          axis.ticks = element_blank(),
          axis.title = element_text(size =15),
          plot.margin=unit(c(15,0,10,0),"mm"),
          legend.text = element_text(size =10),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12, face="bold"))
}







