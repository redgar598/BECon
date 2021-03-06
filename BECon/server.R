library(ggplot2)
library(reshape)


options(shiny.maxRequestSize=30*1024^2)

load("./data/BLBR_app_Objects_beta.RData")
#load("~/BECon/data/BLBR_app_Objects_beta.RData")


# ## density plots stuff
# correlations_BLBR_densityplt<-correlations_BLBR
# correlations_BLBR_densityplt$mean<-rowMeans(correlations_BLBR_densityplt[,2:4])
# 
# var_tissues_densityplt<-var_tissues[,1:3]
# var_tissues_densityplt<-melt(var_tissues_densityplt, id="CpG")
# colnames(var_tissues_densityplt)<-c("CpG","Tissue","Variability")
# levels(var_tissues_densityplt$Tissue)<-c("Blood","All Brain")

#### input summary
CpG_list<-function(gene_interest, CpGs){
  if(is.character(gene_interest)){ 
    Genes<-lapply(1:length(gene_interest), function(x){
      Genes<-Gene_CpG_Relations_minimal[c(grep(paste(", ",gene_interest[x],sep=""), Gene_CpG_Relations_minimal$gene),
                                          grep(paste(gene_interest[x],",",sep=""), Gene_CpG_Relations_minimal$gene),
                                          which(Gene_CpG_Relations_minimal$gene==gene_interest[x]),
                                          which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs)),]
      list_gene<-as.list(Genes$gene)
      list_gene<-sapply(list_gene, function(x) strsplit(x, ", "))
      inlist<-unlist(sapply(1:length(list_gene), function(genes) if(any(list_gene[[genes]]==gene_interest[x])){genes}else{}))
      Genes<-Genes[inlist,]
      Genes})
    Genes<-do.call(rbind, Genes)
    
  }else{
    Genes<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs),]}
  
  Genes_onarray<-unique(Genes$gene)
  Genes_onarray<-unique(unlist(strsplit(Genes_onarray, ",| ,|, ")))
  CpGs<-unique(c(as.character(Genes$Probe_ID), CpGs)) 
  withcorrelation<-nrow(correlations_BLBR[which(correlations_BLBR$CpG%in%CpGs),])
  paste(length(CpGs), "CpGs associated with ", length(Genes_onarray), "genes to look at (",
        withcorrelation, "CpGs have correlation data between blood and brain).", sep=" ")
}


### make CpG list for other inputs
CpG_list_forplottable<-function(gene_interest, CpGs){
  
  if(is.character(gene_interest)){ 
    Genes<-lapply(1:length(gene_interest), function(x){
      Genes<-Gene_CpG_Relations_minimal[c(grep(paste(", ",gene_interest[x],sep=""), Gene_CpG_Relations_minimal$gene),
                                          grep(paste(gene_interest[x],",",sep=""), Gene_CpG_Relations_minimal$gene),
                                          which(Gene_CpG_Relations_minimal$gene==gene_interest[x]),
                                          which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs)),]
      list_gene<-as.list(Genes$gene)
      list_gene<-sapply(list_gene, function(x) strsplit(x, ", "))
      inlist<-unlist(sapply(1:length(list_gene), function(genes) if(any(list_gene[[genes]]==gene_interest[x])){genes}else{}))
      Genes<-Genes[inlist,]
      Genes})
    Genes<-do.call(rbind, Genes)
    
  }else{
    Genes<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs),]}
  CpGs<-unique(c(as.character(Genes$Probe_ID), CpGs)) 
  CpGs
}


################## define drop downs
CpG_names<-as.character(rownames(combat_BLBR_Beta_adjusted))
#gene_names<-as.character(unique(Gene_CpG_Relations_minimal$gene))
gene_names<-Gene_CpG_Relations_update_genes


############################## Comethylation Plot
Comethylation_Plot<-function(correlations_BLBR, Betas, CpG_Hit_List, MaxCpGNum){
  hits_BLBR<-correlations_BLBR[which(correlations_BLBR$CpG%in%CpG_Hit_List),]
  hits_BLBR<-merge(hits_BLBR, Gene_CpG_Relations_minimal[,1:3],  by.x="CpG",by.y="Probe_ID")
  hits_BLBR<-hits_BLBR[with(hits_BLBR, (order(Chromosome_37, Coordinate_37))), ]
  
  if(nrow(hits_BLBR)>MaxCpGNum){hits_BLBR<-hits_BLBR[1:MaxCpGNum,]}else{}
  
  nrow(hits_BLBR)
  hits_BLBR_PBMC_correlation_melt<-melt(hits_BLBR, id="CpG")
  
  meta$SubNumber<-as.factor(meta$SubjectNumC)
  levels(meta$SubNumber)<-c(1:16)
  BLBR_Beta<-as.data.frame(Betas[which(rownames(Betas)%in%hits_BLBR$CpG),])
  BLBR_Beta$CpG<-rownames(BLBR_Beta)
  
  Gene_CpG<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%CpG_Hit_List),c(1:4)]
  Gene_CpG<-Gene_CpG[!duplicated(Gene_CpG),]
  Gene_CpG$Gene_CpG<-sapply(1:nrow(Gene_CpG), function(x) paste(Gene_CpG$Probe_ID[x], " - ",unique(Gene_CpG$gene[x]), sep=""))
  BLBR_Beta<-melt(BLBR_Beta)
  BLBR_Beta<-merge(BLBR_Beta, meta, by.x="variable", by.y="X")
  BLBR_Beta<-merge(BLBR_Beta, Gene_CpG, by.x="CpG", by.y="Probe_ID")
  BLBR_Beta<-BLBR_Beta[with(BLBR_Beta, (order(Chromosome_37, Coordinate_37))),]
  BLBR_Beta$Gene_CpG<-as.factor(BLBR_Beta$Gene_CpG)
  
  BLBR_Beta_order<-BLBR_Beta[,c("Gene_CpG","Coordinate_37","Chromosome_37")]
  BLBR_Beta_order<-BLBR_Beta_order[!duplicated(BLBR_Beta_order), ]
  BLBR_Beta_order<-BLBR_Beta_order$Gene_CpG[with(BLBR_Beta_order, (order(Chromosome_37, Coordinate_37)))]
  
  BLBR_Beta$Gene_CpG<-factor(BLBR_Beta$Gene_CpG, levels=as.character(BLBR_Beta_order))
  levels(BLBR_Beta$TissueType)<-c("Brodmann area 10", "Brodmann area 20", "Brodmann area 7", "Blood")
  
  ggplot()+geom_line(aes(SubNumber, value,group=TissueType,color=TissueType),#,alpha=line
                     BLBR_Beta, size=1.5)+
    theme_bw()+facet_wrap(~Gene_CpG)+
    scale_color_manual(values=c("#fb6a4a","#ef3b2c","#cb181d","cornflowerblue"))+
    ylim(0,1)+xlab("Subject Number")+ylab("Beta Value")+
    theme(text = element_text(size=15))
}

############################### SUmmary Table
summTable<-function(CpGs){
  ## associated Gene information
  CpG_gene<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs),]
  
  ##correlation information
  correlations<-correlations_BLBR[which(correlations_BLBR$CpG%in%CpGs),]
  CpG_gene_correlations<-merge(CpG_gene, correlations, by.x="Probe_ID", by.y="CpG")
  
  ##varbility information
  varibility<-var_tissues[which(var_tissues$CpG%in%CpGs),]
  CpG_gene_correlations_varibility<-merge(CpG_gene_correlations, var_tissues, by.x="Probe_ID", by.y="CpG")
  CpG_gene_correlations_varibility<-CpG_gene_correlations_varibility[order(CpG_gene_correlations_varibility$Probe_ID),]
  
  CpG_gene_correlations_varibility$Coordinate_37<-as.numeric(as.character(CpG_gene_correlations_varibility$Coordinate_37))
  CpG_gene_correlations_varibility<-CpG_gene_correlations_varibility[with(CpG_gene_correlations_varibility, order(gene, Coordinate_37)), ]
  
  CpG_gene_correlations_varibility<-as.data.frame(lapply(CpG_gene_correlations_varibility, as.character))
  CpG_gene_correlations_varibility<-CpG_gene_correlations_varibility[,c(1:12,17:25)]
  
  colnames(CpG_gene_correlations_varibility)<-c("CpG ID","Genomic Coordinate (hg19)","Chr (hg19)",
                                                "Associated Genes", "CpG in Feature of Gene, respectively","Cor Blood-BA7",
                                                "Cor Blood- BA10","Cor Blood- BA20","Mean Cor All Brain",
                                                "SD Mean Cor All Brain", "Percentile of Mean Cor All Brain (positive)",
                                                "Percentile of Mean Cor All Brain (negative)",
                                                "Mean Change Beta with Blood Cell Composition Normalization",
                                                "Mean Change Beta with Brain Cell Composition Normalization",
                                                "Percentile rank of CpG Change Beta with Blood Cell Comp",
                                                "Percentile rank of CpG Change Beta with Brain Cell Comp",
                                                "Var in Blood",
                                                "Var in All Brain","Var in BA7","Var in BA10","Var in BA20")
  
  CpG_gene_correlations_varibility
}

############################### Summary plot (instead of table)


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
fillscale <- scale_fill_manual(name = "Correlation or Cell Composition Percentile\n or Variability Status",values = myColors, drop = FALSE)



Colorful_table_plot<-function(CpGs){
  color_table_plot_cor<-correlations_BLBR[which(correlations_BLBR$CpG%in%CpGs),c(1:4,10:16)]
  color_table_plot_var<-var_tissues[which(var_tissues$CpG%in%CpGs),c(1,2,4:6,8:11)]
  
  color_table_plot_cor_values<-melt(color_table_plot_cor[,c(1:4,8,9)])
  color_table_plot_var_values<-melt(color_table_plot_var[,c(1:5)])
  
  color_table_plot_cor_color<-melt(color_table_plot_cor[,c(1,5:7,10:11)], id="CpG")
  color_table_plot_var_color<-melt(color_table_plot_var[,c(1,6:9)], id="CpG")
  
  color_table_plot_values<-rbind(color_table_plot_cor_values, color_table_plot_var_values)
  color_table_plot_values$Data<-as.factor(color_table_plot_values$variable)
  levels(color_table_plot_values$Data)<-c("correlations","correlations","correlations","Cell_composition","Cell_composition","variability","variability","variability","variability")
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
  levels(cor_var_for_plot$Data)<-c("Correlation","Cell \nComposition","Variability","Chr","Coor","Gene(s)","Gene \nRegion(s)")
  
  cor_var_for_plot$Data<-factor(cor_var_for_plot$Data, levels=c("Chr","Coor","Gene(s)","Gene \nRegion(s)",
                                                                "Variability","Correlation","Cell \nComposition"))
  
  cor_var_for_plot$CpG<-as.factor(cor_var_for_plot$CpG)
  #CpG order
  CpG_order<-Gene_CpG_Relations_minimal[which(Gene_CpG_Relations_minimal$Probe_ID%in%CpGs),]
  CpG_order<-CpG_order[with(CpG_order, rev(order(Chromosome_37, Coordinate_37))), ]
  
  cor_var_for_plot$CpG<-factor(cor_var_for_plot$CpG, 
                               levels=CpG_order$Probe_ID)
  
  # #width
  cor_var_for_plot$variable<-as.character(cor_var_for_plot$variable)
  cor_var_for_plot$variable<-sapply(1:nrow(cor_var_for_plot), function(x) {
    if(cor_var_for_plot$Data[x]=="Chr"){"  "}else{
      if(cor_var_for_plot$Data[x]=="Coor"){"   "}else{
        cor_var_for_plot$variable[x]}}
  })
  
  cor_var_for_plot$variable<-as.factor(cor_var_for_plot$variable)
  cor_var_for_plot$w<-cor_var_for_plot$variable
  levels(cor_var_for_plot$w)<-c(4,0.6,2.8,1,1,1,1,1)
  cor_var_for_plot$w<-as.numeric(as.character(cor_var_for_plot$w))
  
  ggplot(cor_var_for_plot, aes(variable, CpG, fill = color, width=w)) + #
    geom_tile(color = "black",size=0.5) +
    theme_gray(8)+ scale_x_discrete(position = "top")+
    facet_grid(.~Data, scales="free_x",space="free_x")+
    geom_text(aes(label=value),color="black")+ #aes(label = round(value,2)), 
    fillscale+xlab("")+ylab("")+
    theme_minimal()+ guides(fill=guide_legend(ncol=2))+
    theme(axis.text.x = element_text(size =10, color="black",margin=margin(0,0,0,0), face="bold"),
          axis.text.y = element_text(size =12, color="black"),
          axis.ticks = element_blank(),
          axis.title = element_text(size =15),
          plot.margin=unit(c(15,0,10,0),"mm"),
          legend.text = element_text(size =10),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12, face="bold"))
  
}





# ############################### correlation information
# cor_density_plot<-function(CpGs){
#   correlations<-correlations_BLBR_densityplt[which(correlations_BLBR_densityplt$CpG%in%CpGs),]
#   
#   ggplot()+geom_density(aes(mean), correlations_BLBR_densityplt, fill="grey75", color="white")+theme_bw()+
#     geom_vline(aes(xintercept=mean), correlations)+
#     theme(text = element_text(size=15))+xlab("CpG Correlation (Mean Across All Brains)")+
#     ylab("Density")
#   }
# 
# ############################### varibility information
# var_density_plot<-function(CpGs){
#   varibility<-var_tissues_densityplt[which(var_tissues_densityplt$CpG%in%CpGs),]
#   
#   ggplot()+geom_density(aes(Varibility, fill=Tissue), var_tissues_densityplt, color="white")+theme_bw()+
#     geom_vline(aes(xintercept=Varibility), varibility)+facet_wrap(~Tissue, ncol=1)+
#     theme(text = element_text(size=15))+xlab("CpG Varibility (Reference Range)")+
#     ylab("Density")+scale_fill_manual(values=c("cornflowerblue","#fb6a4a"))}
# 












shinyServer(function(input, output, session) {
  
  ## trying to get the CpG list drop down menu 
  
  updateSelectizeInput(session, 'CpG_list', choices = CpG_names, server = TRUE)
  
  ## trying to get the gene list drop down menu 
  
  updateSelectizeInput(session, 'gene_list', choices = gene_names, server = TRUE) 
  
  
  
  
  ## tryin to get CSV upload to work
  filedata <- reactive({
    infile <- input$datafile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.csv(infile$datapath, header=F)
  })
  
  
  
  # helpful text output
  output$text2 <- renderText({
    df <-filedata()
    #if (is.null(df)) return(NULL)
    csvcpg=as.character(df[,1])
    
    paste("You have selected", CpG_list(input$gene_list, c(input$CpG_list, csvcpg))) })
  #   
  

  
  
  ####### COMETHYLATION PLOTS
  output$plot1<-renderPlot({
    df <-filedata()
    #if (is.null(df)) return(NULL)
    csvcpg=as.character(df[,1])
    
    if(length(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))==0){
    plt<-Comethylation_Plot(correlations_BLBR,combat_BLBR_Beta_adjusted, c("cg00308130","cg00201133"), input$CpGnum)
    print(plt)}else{
      plt<-Comethylation_Plot(correlations_BLBR,combat_BLBR_Beta_adjusted, CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)), input$CpGnum)
      print(plt)}
  })
  
 
  
  # Download plots
  plt_dwn<-reactive({
    df <-filedata()
    #if (is.null(df)) return(NULL)
    csvcpg=as.character(df[,1])
    
    if(length(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))==0){
      plt<-Comethylation_Plot(correlations_BLBR,combat_BLBR_Beta_adjusted, c("cg00308130","cg00201133"), input$CpGnum)
      print(plt)}else{
        plt<-Comethylation_Plot(correlations_BLBR,combat_BLBR_Beta_adjusted, CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)), input$CpGnum)
        print(plt)}
  })
  output$downloadPlot <- downloadHandler(
    filename = function() {paste('comethylationplots-', Sys.Date(), '.pdf', sep='')},
    content = function(file) {
      pdf(file, width=18, height=9)
      plt_dwn()
      dev.off()
    })
  
  
  
  
  
  ############## Summary Table Render
  output$view <- renderTable({
    df <-filedata()
    #if (is.null(df)) return(NULL)
    csvcpg=as.character(df[,1])
    
    if(length(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))==0){
      summTable(c("cg00308130","cg00201133"))
    }else{
      summTable(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))}}, include.rownames=FALSE)
  
  # # Download table
  tab<-reactive({
    df <-filedata()
    #if (is.null(df)) return(NULL)
    csvcpg=as.character(df[,1])
    
    if(length(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))==0){
    summTable(c("cg00308130","cg00201133"))
  }else{
    summTable(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))}})
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(tab(), con)
    })
  
  
  #TABLE PLOT
  ## define height of plots depenedent on the CpGs used
  plotHeight = reactive({
    df <-filedata()
    #if (is.null(df)) return(NULL)
    csvcpg=as.character(df[,1])
    
    
    if (length(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg))) < 10) 
      return(300)
    length(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg))) * 40
  })
  
  output$plotsum<-renderPlot({
    
    df <-filedata()
    #if (is.null(df)) return(NULL)
    csvcpg=as.character(df[,1])
    
    if(length(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))==0){
    plt<-Colorful_table_plot(c("cg00308130","cg00201133"))
    print(plt)}else{
      plt<-Colorful_table_plot(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))
      print(plt)}
  }, height=plotHeight,units="px")
  
  
  # Download plots
  plt_tbl_dwn<-reactive({
    df <-filedata()
    #if (is.null(df)) return(NULL)
    csvcpg=as.character(df[,1])
    
    if(length(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))==0){
    plt<-Colorful_table_plot(c("cg00308130","cg00201133"))
    print(plt)}else{
      plt<-Colorful_table_plot(CpG_list_forplottable(input$gene_list,c(input$CpG_list, csvcpg)))
      print(plt)}
  })
  
  
  output$downloadPlotTbl <- downloadHandler(
    filename = function() {paste('SummaryTable-', Sys.Date(), '.pdf', sep='')},
    content = function(file) {
      pdf(file, width=17)
      plt_tbl_dwn()
      dev.off()
    })
  
  
  
  
  #   ######### Correlation density plot
  #   output$plot2<-renderPlot({if(length(CpG_list_forplottable(input$gene_list,input$CpG_list))==0){
  #     plt<-cor_density_plot(c("cg00308130","cg00201133"))
  #     print(plt)}else{
  #       plt<-cor_density_plot(CpG_list_forplottable(input$gene_list,input$CpG_list))
  #       print(plt)}
  #   })
  #   
  #   ########## Variability density plot
  #   output$plot3<-renderPlot({if(length(CpG_list_forplottable(input$gene_list,input$CpG_list))==0){
  #     plt<-var_density_plot(c("cg00308130","cg00201133"))
  #     print(plt)}else{
  #       plt<-var_density_plot(CpG_list_forplottable(input$gene_list,input$CpG_list))
  #       print(plt)}
  #   })
  #   
  
  
  
  ##output$plot1<-renderPlot(plt<-Comethylation_Plot(correlations_BLBR,combat_BLBR_Beta_adjusted, CpG_list_forplottable(input$gene_list,input$CpG_list), input$CpGnum)
  #print(plt)})
  #output$view <- renderTable({
  #summTable(CpG_list_forplottable(input$gene_list,input$CpG_list))}, include.rownames=FALSE)
  
  
  # output$view <- renderTable({
  #summTable(CpG_list_forplottable(input$genename,input$cpgname))}, include.rownames=FALSE)
  
  #output$plot1<-renderPlot({plt<-Comethylation_Plot(correlations_BLBR,combat_BLBR_Beta_adjusted, correlations_BLBR$CpG[1:input$CpGnum])
  #print(plt)})
  
  # Show the first "n" observations
  #output$view <- renderTable({
  #summTable(correlations_BLBR$CpG[1:input$CpGnum])}, include.rownames=FALSE)
  
  
})
