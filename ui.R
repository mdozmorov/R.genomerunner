library(shiny)
library(d3heatmap)

shinyUI( 
  tabsetPanel(
    
    tabPanel("Heatmap",
      fluidPage(
        fluidRow(
          column(4,
                 selectInput("cmbHeatmap", label = "Select which matrix to visualize", 
                             choices = list("P-value" = "matrix_PVAL.txt", 
                                            "Odds Ratio" = "matrix_OR.txt"))),
          column(4,
                 selectInput('cmbCorType',label = "Correlation coefficient type",
                             choices = list("Pearson's" = "pearson",
                                            "Spearman's" = "spearman"))
          ),
          
          column(4,
                 selectInput("cmbClustMethod",label = "Clustering method (hclust)", 
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
        d3heatmapOutput("heatmap", width = "100%", height = "400px"),
        fluidRow(
          column(12,
                 sliderInput("sldNumClust","Number of clusters",min = 2,max=10,value = 3))
        ),
        plotOutput("pltDend",width = "100%", height = "500px")
      )),
    tabPanel("Epigenetic results",
             fluidPage(h3("test"))
             )
    )
  )
