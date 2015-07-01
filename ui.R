library(shiny)
library(d3heatmap)
library(DT)

shinyUI(fluidPage(
    mainPanel(
      tags$head(tags$style("td {align:center;}")),
      uiOutput("mainpage")
    )
))
