library(shiny)
library(d3heatmap)
library(DT)

shinyUI(
  tabsetPanel(
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
             d3heatmapOutput("heatmapEnrich", width = "100%", height = "600px"),
             plotOutput("legendEnrich",width="300px",height="200px")
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
)
