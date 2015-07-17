suppressMessages(source("utils2.R")) # See the required packages there
suppressMessages(source("episimilarity.R"))
(source("genomerunner_d3heatmap.R"))
library(d3heatmap)
library(dendextendRcpp) # required for extracting the height from the dendrogram
library(tools)
library(colorRamps)
library(shinyBS)

# # Lukas paths
results.dir <- "/home/lukas/db_2.00_06-10-2015/results/"
# Mikhail paths
#results.dir <- "/Users/mikhail/Documents/Work/WorkOMRF/Dennis/data.1/chromStates18/"


genomerunner.mode <- TRUE
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
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep=""))
  })
  
  get.adjust.matrix <- reactive({
    mtx <- get.matrix()
    mtx.adjust <- mtx.transform(apply(mtx.untransform(mtx), 2, function(x) p.adjust(abs(x), method = input$cmbPvalAdjustMethod)))
    mtx.adjust <- sign(mtx)*mtx.adjust
  })
  
  output$heatmapEnrich <- renderD3heatmap({
    mtx <- get.adjust.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    d3heatmap::d3heatmap(as.matrix(mtx),heatmap_options = list(hclust=function(tmp) {hclust(tmp, method = input$cmbClustMethod)}), colors = coloring(coloring.num), show_tip=FALSE,dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""))
    
  })
  
  output$legendEnrich <- renderPlot({
    mtx <- get.adjust.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    color.range = seq(min(mtx),max(mtx), (max(mtx)-min(mtx))/coloring.num )
    plot(color.range,rep(1,coloring.num+1),col=coloring(coloring.num+1),pch=15,cex=10,main="Depletion/Enrichment Significance",ylab="",xlab="",yaxt="n")
  })
  
  # generate the enrichment table and appends the gf.name information columns
  get.enrichment.table <- reactive({
    mtx <- read.csv(paste(get.results.dir(),input$cmbMatrix,sep = ""),sep="\t")
    selectedFOI <- 1
    selectedFOI <-input$cmbFOI
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      mtx.adjust <- apply(mtx[selectedFOI], 2, function(x) p.adjust(abs(x), method = input$cmbPvalAdjustMethod))
      pval.sig <- rep("Not significant",nrow(mtx.adjust))
      pval.sig[ mtx[, selectedFOI] < 0.05 & mtx[, selectedFOI] > 0] <- "Overrepresented" 
      pval.sig[mtx[, selectedFOI] > -0.05 & mtx[, selectedFOI] < 0] <- "Underrepresented"
      
      mtx.table <- cbind(abs(mtx[selectedFOI]), 
                         pval.sig, 
                         mtx.adjust)
      colnames(mtx.table) <- c("P.value","Direction","P.adj")
    }else{ 
      # for odds ratio
      or.sig <- rep("Not significant", nrow(mtx))
      or.sig[mtx[, selectedFOI] > 1] <- "Enriched"
      or.sig[mtx[, selectedFOI] < 1] <- "Depleted"
      mtx.table <- cbind(abs(mtx[selectedFOI]),
                         or.sig)
      colnames(mtx.table) <- c("Odds.ratio","Direction")
    }
    mtx.table = cbind(GF.Name = rownames(mtx.table),mtx.table) # make the GF.name a column instead of just the rowname so we can left_join
    rownames(mtx.table) <- NULL
    # gfAnnot is loaded in utils2
    mtx.table <- left_join(mtx.table,gfAnnot,by=c("GF.Name"="name"))
    mtx.table <- subset(mtx.table, select = -c(ind))
    
    return(mtx.table)
  })
  
  output$tblEnrichment <-renderDataTable({
    get.enrichment.table()
  }, options = list( lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                     pageLength = 50))
  
  output$downloadEnrichTable <- downloadHandler(
    filename = function() { 
      return("Enrichment_table.txt")
    },
    content = function(file) {
      write.table(x = get.enrichment.table(),file =  file ,sep = "\t",quote = F,row.names = F)
    }
  )
  
  
  outputOptions(output, "downloadEnrichTable", suspendWhenHidden=FALSE)
  
  ## enrichment up and down plots for single column --------------------------------------------------------------
  get.barplot.matrix <- reactive({
    # populate the enrichment table combobox
    file.names.enrichment <- file_path_sans_ext(list.files(paste(get.results.dir(),"enrichment/",sep="")))
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep=""))
    
  })
  
  
  getEnrichmentUp <- reactive({
    mtx <-  data.frame( get.barplot.matrix())
   
    
    selectedFOI <- 1
    selectedFOI <-input$cmbFOI
    updown.split =  0
    mtx.up <- subset(mtx[selectedFOI],mtx[selectedFOI] > updown.split)
    
    if (nrow(mtx.up)==0){
      return(mtx.up)
    }
    
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      
      # do the pvalue adjustment 
      mtx.up.adjust <- apply(1/10^(abs(mtx.up)), 2, function(x) {p.adjust(x,  method = input$cmbPvalAdjustMethod)})
      mtx.up.adjust <- as.matrix(mtx.up.adjust)
      rownames(mtx.up.adjust) <- rownames(mtx.up); colnames(mtx.up.adjust) <- colnames(mtx.up);
      mtx.up <- mtx.transform(mtx.up.adjust);
    }
    
    
    # sort the results
    mtx.up.sorted <- mtx.up[order(mtx.up,decreasing = T),,drop=FALSE]
  })
  
  getEnrichmentDown <- reactive({
    mtx <-  data.frame( get.barplot.matrix())
    selectedFOI <- 1
    selectedFOI <-input$cmbFOI
    updown.split = 0
    mtx <-  data.frame(get.barplot.matrix())
    mtx.down <- subset(mtx[selectedFOI],mtx[selectedFOI] < updown.split,drop = F)
    
    if (nrow(mtx.down)==0){
      return(mtx.down)
    }
    
    if (input$cmbMatrix == "matrix_PVAL.txt"){
      
      # do the pvalue adjustment 
      mtx.down.adjust <- apply(1/10^(abs(mtx.down)), 2, function(x) {p.adjust(x,  method = input$cmbPvalAdjustMethod)})
      mtx.down.adjust <- -as.matrix(mtx.down.adjust) # apply negative bc everything is downregulated here
      rownames(mtx.down.adjust) <- rownames(mtx.down); colnames(mtx.down.adjust) <- colnames(mtx.down);
      mtx.down <- mtx.transform(mtx.down.adjust);
    }
    mtx.down.sorted <- mtx.down[order(mtx.down,decreasing = F),,drop=FALSE]
  })
  
  # the same barplot is used for single FOI results and multiple FOI results
  output$pltEnrichUp <- renderPlot({
    mtx.up.sorted = getEnrichmentUp()
    if (nrow(mtx.up.sorted)==0){
      # plot raw text
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("Nothing overrepresented to plot."),  cex = 1.6, col = "black")
      box()
      return()
    }
    par(mar = c(10,5,4.1,2.1))
    if (nrow(mtx.up.sorted) < 3){
      barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
              space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.up.sorted),30),main="Enriched epigenomic associations",xlim=c(0,10*nrow(mtx.up.sorted)),axes=F)
    } else{
      barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
              space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.up.sorted),30),main="Enriched epigenomic associations",axes=F)
    }
    abline(a=0,b=0)   
    if(input$cmbMatrix == "matrix_PVAL.txt"){
      mtext('P-value',side=2,line=4)
      axis.values = seq(0,100,by=1)
      axis(side = 2, at = axis.values,labels=1/10^abs(axis.values),las=1)
    } else{
      mtext('Odds-ratio',side=2,line=4)
      axis.values = seq(0,100,by = .1)
      axis(side = 2, at = axis.values, labels=round(2^axis.values,digits = 2),las=1)
    }
  })
  
  output$pltEnrichDown <- renderPlot({
    mtx.down.sorted <- abs(getEnrichmentDown())
    if (nrow(mtx.down.sorted) == 0){
      # plot raw text
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("Nothing underrepresented to plot."),  cex = 1.6, col = "black")
      box()
      return(mtx.down.sorted)
    }
    par(mar = c(10,5,4.1,2.1))
    if (nrow(mtx.down.sorted) < 3){
      barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
              space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),main = "Depleted epigenomic associations",xlim=c(0,10*nrow(mtx.down.sorted)),axes = F)
    }else{
      barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
              space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),main = "Depleted epigenomic associations", axes = F)
    }
    abline(a=0,b=0)   
    if(input$cmbMatrix == "matrix_PVAL.txt"){
      mtext('P-value',side=2,line=4)
      axis.values = c(seq(0,100,by=1))
      axis(side = 2, at = axis.values,labels=1/10^abs(axis.values), las=1)
    }
    else{
      mtext('Odds-ratio',side=2,line=4)
      axis.values = c(seq(0,100, by=.5))
      axis(side = 2, at = axis.values, labels=round(2^axis.values,digits = 2), las=1)
    }
  })
  
  # episimilarity ---------------------------------------------------------------
  get.corr.matrix <- reactive({
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep=""))
    mtx <- scale(mtx)
    
    rcorr(as.matrix(mtx), type=input$cmbEpisimCorType)[[1]]
  })
  
  get.cor.hclust.dendrogram <- reactive({
    cor.mat <- get.corr.matrix() 
    as.dendrogram(hclust(as.dist(1-cor.mat), method=input$cmbClustMethod))
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
    plot(color.range,rep(1,coloring.num+1),col=coloring(coloring.num+1),pch=15,cex=10,main="Correlation Coefficient",ylab="",xlab="",yaxt="n")
  })
  
  # this function is cut out from the tblEpigenetics renderer. It is a long calculation that is only run when # of clusters changes
  calculate.clust <- reactive({
    cor.mat <- get.corr.matrix() # this line ensure that dendrogram is redrawn when heatmap is
    hclustergram <- get.cor.hclust.dendrogram() # ensures that dendrogram is redrawn when hclust method is changed
    
    dend = readRDS(file = paste(get.results.dir(), "heatmap.dend.rds",sep=""))
    cl_num <- input$sldEpisimNumClust # Empirically set desired numter of clusters
    hcut <- heights_per_k.dendrogram(dend)[cl_num] # extract the height to cut based on # of groups
    # get the cluster labels
    mtx.clust <- dend %>% mtx.clusters(height=hcut, minmembers=3)
    mtx = load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep="")) # load the original matrix
    mtx.deg <- suppressWarnings(mtx.degfs(mtx[, mtx.clust$eset.labels], mtx.clust, label="broadPeak2"))
    
    updateSelectInput(session,"cmbEpigenetics","Select which epigenetic table to render",choices = names(mtx.deg))
    mtx.deg.path = paste(get.results.dir(),"mtx.deg.episim.RDS",sep = "")
    saveRDS(mtx.deg,file=mtx.deg.path)
    return(mtx.deg)
  })
  
  get.epigenetics.table <- reactive({
    mtx.deg <- calculate.clust()
    # check if any results were returned
    if (is.null(names(mtx.deg))){ 
      return(data.frame(NoResult="There is nothing signficant to show"))
    }
    selectedCor = names(mtx.deg)[1]
    # save last selected value
    if (input$cmbEpigenetics != "Results not ready yet.") {
      if (input$cmbEpigenetics %in% names(mtx.deg)){
        selectedCor = input$cmbEpigenetics
      }
      updateSelectInput(session,"cmbEpigenetics",choices=names(mtx.deg),selected = selectedCor)
    }
    #convert values to numeric form for sorting purposes
    for(x in list("adj.p.val",3,4)){
      mtx.deg[[selectedCor]][[x]] <- as.numeric(mtx.deg[[selectedCor]][[x]])
    }
    mtx.deg[[selectedCor]] <- subset( mtx.deg[[selectedCor]], select = -c(ind))
    mtx.deg[[selectedCor]]
  })
  
  output$tblEpigenetics <-renderDataTable({
    validate(need(try(get.epigenetics.table()),"Try a different clustering method."))
    get.epigenetics.table()
  },options = list( lengthMenu = list(c(10, 50, 100,-1), c('10', '50','100', 'All')),
                    pageLength = 50))
  
  output$downloadEpigenetics <- downloadHandler(
    filename = function() { 
      return("Epigenetics_table.txt")
    },
    content = function(file) {
      write.table(x = get.epigenetics.table(),file =  file ,sep = "\t",quote = F,row.names = F)
    }
  )
  
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
    mtx = load_gr_data(paste(get.results.dir(), input$cmbMatrix,sep="")) # load the original matrix
    
    # Define the clusters by rectangles
    validate(need(try(rect.hclust( as.hclust(dend), k=cl_num, border=cols)),"Try a different clustering method."))
    
  })
  
  
  # --download button code --------------------------------------------------
  output$downloadEnrichBarPDF <- downloadHandler(
    filename = function() { 
      return("EnrichmentUp.pdf")
    },
    content = function(file) {
      pdf(file=file,width=9,height=5)
      # print EnrichmentUp to PDF
      mtx.up.sorted = getEnrichmentUp()
      if (nrow(mtx.up.sorted)==0){
        # plot raw text
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("Nothing overrepresented to plot."),  cex = 1.6, col = "black")
        box()
      }else{
        par(mar = c(10,5,4.1,2.1))
        if (nrow(mtx.up.sorted) < 3){
          barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.up.sorted),30),main="Enriched epigenomic associations",xlim=c(0,10*nrow(mtx.up.sorted)),axes=F)
        } else{
          barplot(as.matrix(t(head(mtx.up.sorted,30))), beside=T,col = "red3",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.up.sorted),30),main="Enriched epigenomic associations",axes=F)
        }
        abline(a=0,b=0)   
        if(input$cmbMatrix == "matrix_PVAL.txt"){
          mtext('P-value',side=2,line=4)
          axis.values = c(seq(0,100,by=1))
          axis(side = 2, at = axis.values,labels=1/10^abs(axis.values),las=1)
        } else{
          mtext('Odds-ratio',side=2,line=4)
          axis.values = seq(0,100,by = 1)
          axis(side = 2, at = axis.values, labels=round(2^axis.values,digits = 2),las=1)
        }
      }
      # print Enrichmentdown bar plot to PDF
      mtx.down.sorted <- abs(getEnrichmentDown())
      if (nrow(mtx.down.sorted)==0){
        # plot raw text
        par(mar = c(0,0,0,0))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(x = 0.5, y = 0.5, paste("Nothing overrepresented to plot."),  cex = 1.6, col = "black")
        box()
      } else{
        par(mar = c(10,5,4.1,2.1))
        if (nrow(mtx.down.sorted) < 3){
          barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),main = "Depleted epigenomic associations",xlim=c(0,10*nrow(mtx.down.sorted)),axes = F)
        }else{
          barplot(as.matrix(t(head(mtx.down.sorted,30))), beside=T,col = "green4",
                  space=c(0.2,1), cex.names=1, las=2, names.arg=head(rownames(mtx.down.sorted),30),main = "Depleted epigenomic associations", axes = F)
        }
        abline(a=0,b=0) 
        if(input$cmbMatrix == "matrix_PVAL.txt"){
          mtext('P-value',side=2,line=4)
          axis.values = c(seq(-100,0,by=1))
          axis(side = 2, at = axis.values,labels=1/10^abs(axis.values), las=1)
        }
        else{
          mtext('Odds-ratio',side=2,line=4)
          axis.values = c(seq(-10, 0, by=.1))
          axis(side = 2, at = axis.values, labels=round(2^axis.values,digits = 2), las=1)
        }
      }
      
      dev.off()
    },
    contentType = 'application/pdf'
  )
  
  
  # Create a different UI depending if there are multiple GF in the results
  output$mainpage <-renderUI({
    mtx <- load_gr_data(paste(get.results.dir(), "matrix_PVAL.txt",sep="")) # manually load matrix since controls are not loaded yet
    
    single.feature = TRUE;
    if (ncol(mtx)>1 & nrow(mtx)>1){single.feature = FALSE}
    if (single.feature == FALSE){
      tabsetPanel(id="tabsMultiple",
                  tabPanel("Enrichment analysis heatmap",
                           d3heatmapOutput("heatmapEnrich", width = "100%", height = "600px"),
                           plotOutput("legendEnrich",width="300px",height="150px")
                  ), 
                  tabPanel("Enrichment analysis barplot",
                           downloadButton('downloadEnrichBarPDF', 'Download PDF'),
                           plotOutput("pltEnrichUp",width="100%",height = "350px"),
                           plotOutput("pltEnrichDown", width="100%", height= "350px")
                  ),
                  tabPanel("Enrichment analysis tables",
                           br(),br(),
                           downloadButton('downloadEnrichTable', 'Download table in tab-separated format'),
                           br(),br(),
                           DT::dataTableOutput("tblEnrichment")),
                  tabPanel("Epigenetic similarity analysis heatmap",
                           fluidPage(
                             fluidRow(
                               column(6,
                                      d3heatmapOutput("heatmapEpisim", width = "100%", height = "600px"),
                                      plotOutput("legendEpisim",width="300px",height="150px")
                               ),
                               column(6, 
                                      plotOutput("pltDend",width = "100%", height = "500px")
                               )
                             )
                           )),
                  tabPanel("Epigenetic similarity analysis tables",
                           selectInput("cmbEpigenetics", "Select which epigenetic analysis to show", choices = list("Results not ready yet.")),
                           br(),br(),
                           downloadButton('downloadEpigenetics', 'Download table in tab-separated format'),
                           br(),br(),
                           DT::dataTableOutput("tblEpigenetics"))
      )
    } else{ # this UI is created when only a single GF result is returned
      tabsetPanel(id="tabsSingleGF",
                  tabPanel("Enrichment analysis barplot",
                           downloadButton('downloadEnrichBarPDF', 'Download PDF'),
                           plotOutput("pltEnrichUp",width="100%",height = "350px"),
                           plotOutput("pltEnrichDown", width="100%", height= "350px")
                  ),
                  tabPanel("Enrichment analysis tables",
                           br(),br(),
                           downloadButton('downloadEnrichTable', 'Download table in tab-separated format'),
                           br(),br(),
                           DT::dataTableOutput("tblEnrichment"))
      )
    }
  })
  output$sidebar <- renderUI({
    mtx <- load_gr_data(paste(get.results.dir(), "matrix_PVAL.txt",sep="")) # manually load matrix since controls are not loaded yet
    
    single.feature = TRUE
    if (ncol(mtx)>1 & nrow(mtx)>1){single.feature = FALSE}
    if (single.feature == FALSE){
      sidebarPanel(width = 4,h3("Global Settings"),
                   selectInput("cmbMatrix", label = "Results to visualize", 
                               choices = list("P-values" = "matrix_PVAL.txt", 
                                              "Odds Ratios" = "matrix_OR.txt")),
                   bsTooltip("cmbMatrix", "Select P-value or Odds ratio", placement = "top", trigger = "hover"),
                   selectInput("cmbFOI", "Select which SNP set to visualize", choices = colnames(mtx)),
                   conditionalPanel("input.cmbMatrix=='matrix_PVAL.txt'",
                                    selectInput("cmbPvalAdjustMethod",label = "P-value multiple testing correction method",
                                                choices = c( "fdr","none","BH","holm", "hochberg", "hommel", "bonferroni","BY"))),
                   conditionalPanel("input.tabsMultiple == 'Enrichment analysis heatmap' || input.tabsMultiple == 'Epigenetic similarity analysis heatmap'",
                                    hr(),h3("Heatmap Settings"),
                                    selectInput("cmbClustMethod",label = "Clustering method (hclust)", 
                                                choices = list("ward.D",
                                                               "ward.D2",
                                                               "single",
                                                               "complete",
                                                               "average",
                                                               "mcquitty",
                                                               "median",
                                                               "centroid")
                                    ),
                                    selectInput('cmbEpisimCorType',label = "Correlation coefficient type",
                                                choices = list("Pearson's" = "pearson",
                                                               "Spearman's" = "spearman"))
                   ),
                   conditionalPanel("input.tabsMultiple == 'Epigenetic similarity analysis heatmap'",
                                    hr(),h3("Epigenetic similarity"),
                                    sliderInput("sldEpisimNumClust","Number of clusters",min = 2,max=10,value = 3)
                   )
                   
      )
    }else{ # this is for a single column result file
      sidebarPanel(h3("Global Settings"), hr(),
                   selectInput("cmbMatrix", label = "Results to visualize", 
                               choices = list("P-values" = "matrix_PVAL.txt", 
                                              "Odds Ratios" = "matrix_OR.txt")),
                   selectInput("cmbFOI", "Select which SNP set to visualize", choices = colnames(mtx)),
                   conditionalPanel("input.cmbMatrix=='matrix_PVAL.txt'",
                                    if(nrow(mtx)>1){
                                      selectInput("cmbPvalAdjustMethod",label = "P-value multiple testing correction method",
                                                  choices = c( "fdr","none","BH","holm", "hochberg", "hommel", "bonferroni","BY"))}
                   )
      )
    }
  })
})

