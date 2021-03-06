shinyUI(fluidPage(
  titlePanel(h1("BECon: A tool for interpreting DNA methylation findings from blood in the context of brain"),
             tags$head(tags$link(rel = "icon", type = "image/png", href = "overview_icon.png"),
                       tags$title("BECon"))),
  
  
  
  
  
  sidebarLayout(
    sidebarPanel( h3("Gene and CpG Selection"),
                  h5("Input CpGs or genes of interest to explore the correlation level between methylation in human blood and brain. 
                     The aim of BECon (Blood-Brain Epigenetic Concordance) is to allow for improved interpretation of surrogate 
                     methylation results by looking at the relationship human blood and brain methylation."),
                  helpText(a("Click Here For More Information",href="Blood_brain_Shiny_Help.pdf", target="_blank")),
                  
                  # Gene Name input
                  selectizeInput('gene_list', h4("Genes To View"), choices=NULL, multiple = TRUE, 
                                 options = list(placeholder = 'Input Gene Names')),
                  
                  # CpG List Input
                  selectizeInput('CpG_list', h4("CpGs To View"), choices=NULL, multiple = TRUE, 
                                 options = list(placeholder = 'Input CpG IDs')),
                  #CpG csv upload
                  fileInput('datafile', h6('Upload CpG IDs From a File \n (see help document for file format)'),
                            accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                  
                  
                  # maximum CpGs to view
                  numericInput("CpGnum", label = h4("Max CpG Number to Plot"), value = 12),
                  
                  submitButton("Update View"),
                  img(src = "overview2.png", height = 320, width = 320),
                  width = 3),
    
    
    
    mainPanel(
      # helpful bit of text
      textOutput("text2"),
      
      
      uiOutput('sized_plot'),
      
      #co meth plot
      downloadLink('downloadPlot', 'Download Plots'),
      plotOutput("plot1", height=400, width = 1000),
      
      
      #       #correlation plot
      #       h4("Correlation Density Plot with Selected CpGs Indicated"),
      #       plotOutput("plot2", height= 400,width=400),
      #       #Varibility plot
      #       h4("Varibility Density Plots with Selected CpGs Indicated"),
      #       plotOutput("plot3", width=400),
      
      # summary table and plot
      h4("Summary of the Correlations and Variability of Selected CpGs"),
      "Download as a", 
      downloadLink('downloadPlotTbl', 'colorful table'), "or a",
      downloadLink('downloadData', 'csv file'),
      plotOutput("plotsum",width = 1200))
    
    #tableOutput("view"))
    # Summary plot
    
    
    
    
    )))
