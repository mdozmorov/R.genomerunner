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
  
  get.hclust.mat <- reactive({
    cor.mat <- get.corr.matrix() 
    hclust(as.dist(1-cor.mat), method=input$cmbClustMethod);
  })
  
  get.dend.file <- reactive({
    readRDS(file = "/home/lukas/heatmap.dend.rds")
  })
  
  output$heatmap <- renderD3heatmap({     
    cor.mat <- get.corr.matrix()
    hclustergram <- get.hclust.mat()
    d3heatmap::d3heatmap(as.matrix(cor.mat),heatmap_options = list(Rowv=hclustergram$Rowv,Colv=hclustergram$Colv,keep.dendro=TRUE))
  })
  output$pltDend <- renderPlot({
    dend = get.dend.file()
    plot(dend)
  })
})
