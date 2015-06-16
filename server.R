suppressMessages(source("utils2.R")) # See the required packages there
suppressMessages(source("episimilarity.R"))
suppressMessages(source("genomerunner_d3heatmap.R"))
library(d3heatmap)
library(dendextendRcpp) # required for extracting the height from the dendrogram

results.dir <- "/home/lukas/db_2.00_06-10-2015/results/test1/"
genomerunner.mode <- FALSE
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
  
  get.corr.matrix <- reactive({
    print(get.results.dir())
    mtx <- load_gr_data(paste(get.results.dir(), input$cmbHeatmap,sep=""))
    mtx <- scale(mtx)
    rcorr(as.matrix(mtx), type=input$cmbCorType)[[1]]
  })
  
  get.hclust.dendrogram <- reactive({
    cor.mat <- get.corr.matrix() 
    as.dendrogram(hclust(as.dist(1-cor.mat), method=input$cmbClustMethod))
  })
  
  output$heatmap <- renderD3heatmap({     
    cor.mat <- get.corr.matrix()
    hclustergram <- get.hclust.dendrogram()
    d3heatmap::d3heatmap(as.matrix(cor.mat),heatmap_options = list(Rowv=hclustergram,Colv=hclustergram,keep.dendro=TRUE))
  })
  
  # Parse the GET query string
  output$queryText <- reactive({
    
  })
  
  output$pltDend <- renderPlot({ 
    cor.mat <- get.corr.matrix() # this line ensure that dendrogram is redrawn when heatmap is
    hclustergram <- get.hclust.dendrogram() # ensures that dendrogram is redrawn when hclust method is changed
    
    dend = readRDS(file = "/home/lukas/heatmap.dend.rds")
    plot(as.dendrogram(dend, hang=-1)) # Plot dendrogram
    cl_num <- input$sldNumClust # Empirically set desired numter of clusters
    cols <- rainbow(cl_num) # Make colors for clusters
    print(cl_num)
    hcut <- heights_per_k.dendrogram(dend)[cl_num] # extract the height to cut based on # of groups
    # get the cluster labels
    mtx.clust <- dend %>% mtx.clusters(height=hcut, minmembers=3)
    write.table(as.data.frame(mtx.clust), "/home/lukas/clustering_all.txt", sep="\t", row.names=FALSE, quote=FALSE)
    mtx = load_gr_data(paste(get.results.dir(), input$cmbHeatmap,sep="")) # load the original matrix
    mtx.deg <- mtx.degfs(mtx[, mtx.clust$eset.labels], mtx.clust, label="broadPeak2")
    rect.hclust(as.hclust(dend), k=cl_num, border=cols) # Define the clusters by rectangles
    # what's the purpose of this code?
    #abline(h=hcut)
    #as.hclust(dend)$height %>% density %>% plot
    #abline(v=hcut)
    
    #mtx.clust <- dend %>% mtx.clusters(height=hcut, minmembers=3)
  })
})
