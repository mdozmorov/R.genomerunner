suppressMessages(source("utils2.R")) # See the required packages there
suppressMessages(source("episimilarity.R"))
(source("genomerunner_d3heatmap.R"))
library(d3heatmap)
library(dendextendRcpp) # required for extracting the height from the dendrogram
library(tools)
library(colorRamps)

# # Lukas paths
# results.dir <- "/home/lukas/db_2.00_06-10-2015/results/test2_single_col/"
# Mikhail paths
results.dir <- "/Users/mikhail/Documents/Work/WorkOMRF/Dennis/data.1/DNAse_hotspotbroadall/"


genomerunner.mode <- FALSE
coloring.num = 50
shinyServer(function(input, output,session) {
  
  
  # Parse the GET query string
  get.results.dir <- reactive({
    if (genomerunner.mode){
      query <- parseQueryString(session$clientData$url_search)
      return(paste(results.dir, query$id,"/",sep = ""))
    }else{
      return(results.dir)
    }
  })
  
  # enrichment ------------------------------------------------------------------
  get.matrix <- reactive({
    # populate the enrichment table combobox
    file.names.enrichment <- file_path_sans_ext(list.files(paste(get.results.dir(),"enrichment/",sep="")))
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbEnrichHeatmap,sep=""))
  })
  
  output$heatmapEnrich <- renderD3heatmap({
    mat <- get.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
  
    if(!check.single_gf()){
      d3heatmap::d3heatmap(as.matrix(mat),heatmap_options = list(hclust=function(tmp) {hclust(tmp, method = input$cmbEnrichClust)}), colors = coloring(coloring.num), show_tip=FALSE,dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""))
    }else{
      NULL
    }
  })
  
  output$legendEnrich <- renderPlot({
    mtx <- get.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    color.range = seq(min(mtx),max(mtx), (max(mtx)-min(mtx))/coloring.num )
    plot(color.range,rep(1,coloring.num+1),col=coloring(coloring.num+1),pch=15,cex=10,main="Heatmap Legend",ylab="",xlab="")
  })
  
  output$tblEnrichment <-renderDataTable({      
    
    mtx <- read.csv(paste(get.results.dir(),input$cmbEnrichBarplot,sep = ""),sep="\t")
    if (input$cmbEnrichBarplot == "matrix_PVAL.txt"){
      mtx.adjust <- apply(mtx, 2, function(x) p.adjust(abs(x), method = input$cmbEnrichBarPlotPvalAdjust))
      mtx.sign <- ifelse(sign(mtx) < 0, "Underrepresented", "Overrepresented")
      
      mtx.table <- cbind(abs(mtx), 
                         mtx.sign, 
                         mtx.adjust)
      colnames(mtx.table) <- c("P.value","Direction","P.adj")
      return(mtx.table)
    }else{ 
      # for odds ratio
      mtx.sign <- ifelse(sign(mtx) < 1, "Underrepresented", "Overrepresented")
      mtx.table <- cbind(abs(mtx),
                         mtx.sign)
      colnames(mtx.table) <- c("Odds.ratio","Direction")
      return(mtx.table)
    }
  })
  
  ## enrichment up and down plots for single column --------------------------------------------------------------
  get.barplot.matrix <- reactive({
    # populate the enrichment table combobox
    file.names.enrichment <- file_path_sans_ext(list.files(paste(get.results.dir(),"enrichment/",sep="")))
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbEnrichBarplot,sep=""))
  })
  
  check.single_gf <- reactive({
    # Check if there is only one column in the matrix.  If so, we will plot a bar plots intead of heatmap
    mtx <- get.barplot.matrix()
    if(ncol(mtx)==1){
      def.value = 30
      if (nrow(mtx)<30){def.value = nrow(mtx)}
      updateSliderInput(session,"sldNumFeatures",min = 1,max = nrow(mtx),value = def.value)
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  
  
  output$pltEnrichUp <- renderPlot({
    if (check.single_gf() != TRUE){
      return(plot.new())
    }
    updown.split =  0
    mtx <-  data.frame( get.barplot.matrix())
    mtx.up <- subset(mtx,mtx[1] > updown.split)
    # filter out results that do not meet pvalue threshold
    log10.pval = -log10(input$numBarplotThreshold)
    if (input$cmbEnrichBarplot == "matrix_PVAL.txt"){
      mtx.up <- subset(mtx.up, mtx.up[1] > log10.pval, drop=F)
      }
    if (nrow(mtx.up)==0){
      # plot raw text
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("Nothing overrepresented to plot."),  cex = 1.6, col = "black")
      box()
      return()
    }
    if (input$cmbEnrichBarplot == "matrix_PVAL.txt"){
      # do the pvalue adjustment 
     #  mtx.up <- p.adjust(mtx.up,method = input$cmbEnrichBarPlotPvalAdjust)
    }
    
    if (input$cmbEnrichBarplot == "matrix_PVAL.txt"){
      
      # do the pvalue adjustment 
      mtx.up.adjust <- apply(1/10^(abs(mtx.up)), 2, function(x) {p.adjust(x,  method = input$cmbEnrichBarPlotPvalAdjust)})
      mtx.up.adjust <- as.matrix(mtx.up.adjust)
      rownames(mtx.up.adjust) <- rownames(mtx.up); colnames(mtx.up.adjust) <- colnames(mtx.up);
      mtx.up <- mtx.transform(mtx.up.adjust);
    }
    
    
    # sort the results
    mtx.up.sorted <- mtx.up[order(mtx.up,decreasing = T),,drop=FALSE]
    if(input$cmbEnrichBarplot == "matrix_PVAL.txt"){
      y.label = "-log10(p-value)"
    }else{y.label="log2(odds-ratio)"}
    par(mar = c(5,5,4.1,2.1))
    barplot(as.matrix(t(head(mtx.up.sorted,input$sldNumFeatures))), beside=T,col = "red3",
            space=c(0.2,1), cex.names=0.8, las=2, names.arg=head(rownames(mtx.up.sorted),input$sldNumFeatures),ylab=y.label,main="Enriched epigenomic associations")
    abline(a=0,b=0)
    
    #barplot1(head(mtx.up.sorted,input$sldNumFeatures),names.args = head(rownames(mtx.up.sorted),input$sldNumFeatures))
    #names.args.up[names.args.up == "NA:NA"] <- make.names(unlist(lapply(mtx.sorted.up, rownames)), unique=T)[names.args.up == "NA:NA"]
    #names.args.dn[names.args.dn == "NA:NA"] <- make.names(unlist(lapply(mtx.sorted.dn, rownames)), unique=T)[names.args.dn == "NA:NA"]
    #bottom <- 8
    # Plot barplots
    #if (!grepl("dn", toPlot) & (nrow(mtx.barplot.up) > 0)) { barplot1(mtx.barplot.up[, seq(1:length(colnum)), drop=F], "topright", bottom=bottom, names.args=names.args.up, pval=pval) }
  })
  
  output$pltEnrichDown <- renderPlot({
    if (check.single_gf() != TRUE){
      return(plot.new())
    }
    
    updown.split = 0
    mtx <-  data.frame(get.barplot.matrix())
    mtx.down <- subset(mtx,mtx[1] < updown.split)
    
    # filter out results that do not meet pvalue threshold
    log10.pval = -log10(input$numBarplotThreshold)
    if (input$cmbEnrichBarplot == "matrix_PVAL.txt"){
      mtx.down <- subset(mtx.down, mtx.down[1] < -log10.pval, drop=F)
    }
    if (nrow(mtx.down)==0){
      # plot raw text
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("Nothing underrepresented to plot."),  cex = 1.6, col = "black")
      box()
      return()
    }   
    
    if (input$cmbEnrichBarplot == "matrix_PVAL.txt"){

      # do the pvalue adjustment 
      mtx.down.adjust <- apply(1/10^(abs(mtx.down)), 2, function(x) {p.adjust(x,  method = input$cmbEnrichBarPlotPvalAdjust)})
      mtx.down.adjust <- -as.matrix(mtx.down.adjust) # apply negative bc everything is downregulated here
      rownames(mtx.down.adjust) <- rownames(mtx.down); colnames(mtx.down.adjust) <- colnames(mtx.down);
      mtx.down <- mtx.transform(mtx.down.adjust);
    }
    mtx.down.sorted <- mtx.down[order(mtx.down,decreasing = F),,drop=FALSE]
    if(input$cmbEnrichBarplot == "matrix_PVAL.txt"){
      y.label = "-log10(p-value)\nnegative = underrepresentation"
    }else{y.label="log2(odds-ratio)\nnegative = underrepresentation"}
    par(mar = c(5,5,4.1,2.1))
    barplot(as.matrix(t(head(mtx.down.sorted,input$sldNumFeatures))), beside=T,col = "green4",
            space=c(0.2,1), cex.names=0.8, las=2, names.arg=head(rownames(mtx.down.sorted),input$sldNumFeatures),ylab=y.label,main = "Depleted epigenomic associations")
    abline(a=0,b=0)
    #barplot(head(mtx.down.sorted,input$sldNumFeatures),names.args = head(rownames(mtx.down.sorted),input$sldNumFeatures))
  })
  
  # episimilarity ---------------------------------------------------------------
  get.corr.matrix <- reactive({
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbEpisimHeatmap,sep=""))
    mtx <- scale(mtx)
    
    rcorr(as.matrix(mtx), type=input$cmbEpisimCorType)[[1]]
  })
  
  get.cor.hclust.dendrogram <- reactive({
    cor.mat <- get.corr.matrix() 
    as.dendrogram(hclust(as.dist(1-cor.mat), method=input$cmbEpisimClustMethod))
  })
  
  output$heatmapEpisim <- renderD3heatmap({     
    cor.mat <- get.corr.matrix()
    hclustergram <- get.cor.hclust.dendrogram()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    d3heatmap::d3heatmap(as.matrix(cor.mat),heatmap_options = list(Rowv=hclustergram,Colv=hclustergram,keep.dendro=TRUE),colors = coloring(coloring.num),dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""))
  })
  
  output$legendEpisim <- renderPlot({
    mtx <- get.corr.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    color.range = seq(min(mtx),max(mtx), (max(mtx)-min(mtx))/coloring.num )
    plot(color.range,rep(1,coloring.num+1),col=coloring(coloring.num+1),pch=15,cex=10,main="Heatmap Legend",ylab="",xlab="")
  })
  
  output$tblEpigenetics <-renderDataTable({
    mtx.deg <- readRDS(file=paste(get.results.dir(),"mtx.deg.episim.RDS",sep = "/"))
    
    #convert values to numeric form for sorting purposes
    for(x in list("adj.p.val",3,4)){
      mtx.deg[[input$cmbEpisimTable]][[x]] <- as.numeric(mtx.deg[[input$cmbEpisimTable]][[x]])
    }
    mtx.deg[[input$cmbEpisimTable]]
  })
  
  output$pltDend <- renderPlot({ 
    cor.mat <- get.corr.matrix() # this line ensure that dendrogram is redrawn when heatmap is
    hclustergram <- get.cor.hclust.dendrogram() # ensures that dendrogram is redrawn when hclust method is changed
    
    dend = readRDS(file = paste(get.results.dir(), "heatmap.dend.rds",sep=""))
    plot(as.dendrogram(dend, hang=-1)) # Plot dendrogram
    cl_num <- input$sldEpisimNumClust # Empirically set desired numter of clusters
    cols <- rainbow(cl_num) # Make colos for clusters
    print(cl_num)
    hcut <- heights_per_k.dendrogram(dend)[cl_num] # extract the height to cut based on # of groups
    # get the cluster labels
    mtx.clust <- dend %>% mtx.clusters(height=hcut, minmembers=3)
    # write.table(as.data.frame(mtx.clust), "/home/lukas/clustering_all.txt", sep="\t", row.names=FALSE, quote=FALSE)
    mtx = load_gr_data(paste(get.results.dir(), input$cmbEpisimHeatmap,sep="")) # load the original matrix
    mtx.deg <- mtx.degfs(mtx[, mtx.clust$eset.labels], mtx.clust, label="broadPeak2")
    
    updateSelectInput(session,"cmbEpigenetics","Select which epigenetic table to render",choices = names(mtx.deg))
    mtx.deg.path = paste(get.results.dir(),"mtx.deg.episim.RDS",sep = "/")
    saveRDS(mtx.deg,file=mtx.deg.path)
    # create the tabs
    
    rect.hclust(as.hclust(dend), k=cl_num, border=cols) # Define the clusters by rectangles
    
  })
  
  # Create a different UI depending if there are multiple GF in the results
  output$mainpage <-renderUI({
    mtx <- load_gr_data(paste(get.results.dir(), "matrix_PVAL.txt",sep="")) # manually load matrix since controls are not loaded yet
    single.gf = TRUE
    if (ncol(mtx)>1){single.gf = FALSE}
    if (single.gf == FALSE){
      tabsetPanel("tabsMultiple",
                  tabPanel("Enrichment analysis barplot",
                           fluidRow(
                             column(6,
                                    selectInput("cmbEnrichBarplot", label = "Select which matrix to visualize", 
                                                choices = list("P-value" = "matrix_PVAL.txt", 
                                                               "Odds Ratio" = "matrix_OR.txt")),
                                    conditionalPanel("input.cmbEnrichBarplot=='matrix_PVAL.txt'",
                                                     selectInput("cmbEnrichBarPlotPvalAdjust",label = "Select which P-value comparison correction",
                                                                 choices = c( "fdr","none","BH","holm", "hochberg", "hommel", "bonferroni","BY")))
                             ),
                             column(6,
                                    conditionalPanel("input.cmbEnrichBarplot=='matrix_OR.txt'",
                                                     sliderInput("sldNumFeatures",label = "Number of top results to plot",min=1,max=100,value=30)),
                                    conditionalPanel("input.cmbEnrichBarplot=='matrix_PVAL.txt'",
                                                     numericInput("numBarplotThreshold","Filter by threshold: lower limit",min = 0,max=1,value = 0.05)))
                           ),
                           plotOutput("pltEnrichUp",width="100%",height = "350px"),
                           plotOutput("pltEnrichDown", width="100%", height= "350px")
                  ),
                  tabPanel("Enrichment analysis heatmap",
                           fluidRow(
                             column(4,
                                    numericInput("numEnrichFilterLower","Filter by threshold: lower limit",min = 2,max=10,value = 3),
                                    numericInput("numEnrichFilterUpper","Filter by threshold: upper limit",min = 2,max=10,value = 3)),
                             column(4,
                                    selectInput("cmbEnrichHeatmap", label = "Select which matrix to visualize", 
                                                choices = list("P-value" = "matrix_PVAL.txt", 
                                                               "Odds Ratio" = "matrix_OR.txt"))),
                             column(4,
                                    conditionalPanel(condition="check.single_gf == true",
                                                     sliderInput("sldNumFeatures",label = "Number of top results to plot",min=1,max=1000,value=30)),
                                    selectInput("cmbEnrichClust",label = "Clustering method (hclust)", 
                                                choices = list("ward.D",
                                                               "ward.D2",
                                                               "single",
                                                               "complete",
                                                               "average",
                                                               "mcquitty",
                                                               "median",
                                                               "centroid")
                                    )
                             )),
                           
                           conditionalPanel(condition="check.single_gf == false",
                                            d3heatmapOutput("heatmapEnrich", width = "100%", height = "600px"),
                                            plotOutput("legendEnrich",width="300px",height="200px")
                           )
                           
                  ), 
                  tabPanel("Enrichment analysis tables",
                           selectInput("cmbEnrichTable","Select which enrichment table to render",choices=list("Enrichment results not ready")),
                           DT::dataTableOutput("tblEnrichment")),
                  tabPanel("Epigenetic similarity analysis heatmap",
                           fluidPage(
                             fluidRow(
                               column(4,
                                      selectInput("cmbEpisimHeatmap", label = "Select which matrix to visualize", 
                                                  choices = list("P-value" = "matrix_PVAL.txt", 
                                                                 "Odds Ratio" = "matrix_OR.txt"))),
                               column(4,
                                      selectInput('cmbEpisimCorType',label = "Correlation coefficient type",
                                                  choices = list("Pearson's" = "pearson",
                                                                 "Spearman's" = "spearman"))
                               ),
                               
                               column(4,
                                      selectInput("cmbEpisimClustMethod",label = "Clustering method (hclust)", 
                                                  choices = list("ward.D",
                                                                 "ward.D2",
                                                                 "single",
                                                                 "complete",
                                                                 "average",
                                                                 "mcquitty",
                                                                 "median",
                                                                 "centroid")
                                      )
                               )),
                             d3heatmapOutput("heatmapEpisim", width = "100%", height = "600px"),
                             plotOutput("legendEpisim",width="300px",height="200px"),
                             fluidRow(
                               column(12,
                                      sliderInput("sldEpisimNumClust","Number of clusters",min = 2,max=10,value = 3))
                             ),
                             plotOutput("pltDend",width = "100%", height = "500px")
                           )),
                  tabPanel("Epigenetic similarity analysis tables",
                           selectInput("cmbEpisimTable","Select which epigenetic table to render",choices=list("Epigenetic results not ready")),
                           DT::dataTableOutput("tblEpigenetics"))
      )
    } else{ # this UI is created when only a single GF result is returned
      tabsetPanel("tabsSingleGF",
                  tabPanel("Enrichment analysis barplot",
                           fluidRow(
                             column(6,
                                    selectInput("cmbEnrichBarplot", label = "Select which matrix to visualize", 
                                                choices = list("P-value" = "matrix_PVAL.txt", 
                                                               "Odds Ratio" = "matrix_OR.txt")),
                                    conditionalPanel("input.cmbEnrichBarplot=='matrix_PVAL.txt'",
                                                     selectInput("cmbEnrichBarPlotPvalAdjust",label = "Select which P-value comparison correction",
                                                                 choices = c( "fdr","none","BH","holm", "hochberg", "hommel", "bonferroni","BY")))
                             ),
                             column(6,
                                    conditionalPanel("input.cmbEnrichBarplot=='matrix_OR.txt'",
                                                     sliderInput("sldNumFeatures",label = "Number of top results to plot",min=1,max=100,value=30)),
                                    conditionalPanel("input.cmbEnrichBarplot=='matrix_PVAL.txt'",
                                                     numericInput("numBarplotThreshold","Filter by threshold: lower limit",min = 0,max=1,value = 0.05)))
                           ),
                           plotOutput("pltEnrichUp",width="100%",height = "350px"),
                           plotOutput("pltEnrichDown", width="100%", height= "350px")
                  ),
                  tabPanel("Enrichment analysis tables",
                           selectInput("cmbEnrichTable","Select which enrichment table to render",choices=list("Enrichment results not ready")),
                           DT::dataTableOutput("tblEnrichment"))
      )
    }
  })
})
