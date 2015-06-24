suppressMessages(source("utils2.R")) # See the required packages there
suppressMessages(source("episimilarity.R"))
(source("genomerunner_d3heatmap.R"))
library(d3heatmap)
library(dendextendRcpp) # required for extracting the height from the dendrogram
library(tools)
library(colorRamps)

#results.dir <- "/Users/mikhail/Documents/Work/GenomeRunner/gwas2bed/autoimmune/R.GR.autoimmune/data.gr/2Roadmap_DNase_narrowPeak/"
results.dir <- "/Users/mikhail/Documents/Work/WorkOMRF/Dennis/data.1/DNAse_hotspotbroadall/"
# results.dir <- "/Users/mikhail/Documents/Work/WorkOMRF/Dennis/data.1/encTfbs/"
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
    # REMOVE 
    #updateSelectInput(session,"cmbEnrichTable","Select which epigenetic table to render",choices = file.names.enrichment)
    mtx <- read.csv(paste(get.results.dir(), input$cmbEnrichHeatmap,sep=""),sep="\t")
    #!!! do we scale?
    # mtx <- scale(mtx) 
  })
  
  output$heatmapEnrich <- renderD3heatmap({
    mat <- get.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    updateNumericInput(session = session,"numEnrichFilterUpper",label=paste("Filter by threshold: Upper value (max=",max(mat),")"),min = min(mat),max = max(mat),value = max(mat))
    updateNumericInput(session = session,"numEnrichFilterLower",label=paste("Filter by threshold: Lower value (min=",min(mat),")"),min = min(mat),max = max(mat),value = min(mat))
    d3heatmap::d3heatmap(as.matrix(mat),heatmap_options = list(hclust=function(tmp) {hclust(tmp, method = input$cmbEnrichClust)}), colors = coloring(coloring.num), show_tip=FALSE,dendro.rds.path=paste(get.results.dir(),"heatmap.dend.rds", sep=""))
  })
  
  output$legendEnrich <- renderPlot({
    mtx <- get.matrix()
    coloring<-colorRampPalette(c("blue", "yellow", "red"))
    color.range = seq(min(mtx),max(mtx), (max(mtx)-min(mtx))/coloring.num )
    plot(color.range,rep(1,coloring.num+1),col=coloring(coloring.num+1),pch=15,cex=10,main="Heatmap Legend",ylab="",xlab="")
  })
  
  output$tblEnrichment <-renderDataTable({      
    
    enrichment.data <- read.csv(paste(get.results.dir(),"enrichment/",input$cmbEnrichTable,".txt",sep = ""),sep="\t")
    #convert values to numeric form for sorting purposes
    for(x in list(2,4)){
      enrichment.data[[x]] <- as.numeric(enrichment.data[[x]])
    }
    enrichment.data
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
})
