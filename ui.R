library(shiny)
library(d3heatmap)

shinyUI(fluidPage(
  titlePanel("GRsnp HeatMap"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("cmbHeatmap", label = "Select which matrix to visualize", 
                  choices = list("P-value" = "/home/lukas/gftest/matrix_PVAL.txt", 
                                 "Odds Ratio" = "/home/lukas/gftest/matrix_OR.txt")),
      selectInput('cmbCorType',label = "Correlation coefficient type",
                  choices = list("Pearson's" = "pearson",
                                 "Spearman's" = "spearman")),
      selectInput("cmbClustMethod",label = "Clustering method (hclust)", 
                  choices = list("ward.D",
                                 "ward.D2",
                                 "single",
                                 "complete",
                                 "average",
                                 "mcquitty",
                                 "median",
                                 "centroid") ),
      sliderInput("sldNumGroups","Number of clusters",min = 0,max=10,value = 3)
      
    ),
    
    mainPanel(d3heatmapOutput("heatmap", width = "100%", height = "400px"),
              plotOutput("pltDend"))
  )
))