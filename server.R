suppressMessages(source("utils2.R")) # See the required packages there
suppressMessages(source("episimilarity.R"))
suppressMessages(source("genomerunner_d3heatmap.R"))
library(d3heatmap)



shinyServer(function(input, output) {
  
  get.corr.matrix <- reactive({
    mtx <- load_gr_data(input$cmbHeatmap)
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
  
  output$pltDend <- renderPlot({ 
    cor.mat <- get.corr.matrix() # this line ensure that dendrogram is redrawn when heatmap is
    hclustergram <- get.hclust.dendrogram() # ensures that dendrogram is redrawn when hclust method is changed
    
    dend = readRDS(file = "/home/lukas/heatmap.dend.rds")
    plot(as.dendrogram(dend, hang=-1)) # Plot dendrogram
    cl_num <- input$sldNumClust # Empirically set desired numter of clusters
    cols <- rainbow(cl_num) # Make colors for clusters
    rect.hclust(as.hclust(dend), k=cl_num, border=cols) # Define the clusters by rectangles
    
    # what's the purpose of this code?
    #abline(h=hcut)
    #as.hclust(dend)$height %>% density %>% plot
    #abline(v=hcut)
    
    #mtx.clust <- dend %>% mtx.clusters(height=hcut, minmembers=3)
  })
})
