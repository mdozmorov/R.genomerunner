library(d3heatmap)
library(dendextendRcpp) # required for extracting the height from the dendrogram
library(tools)
library(colorRamps)
library(Hmisc)
library(shinyBS)
library(scales)
library(dplyr)
library(gplots)

gfAnnot <- read.table("data/gf_descriptions.txt",sep="\t",header=T, stringsAsFactors = FALSE)
