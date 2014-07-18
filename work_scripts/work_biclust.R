

# === Biclustering
library(biclust)
library(eisa)
library("ExpressionView")
library(Biobase)
mtx.eset<-new("ExpressionSet",exprs=as.matrix(mtx))
Bc <- biclust(as.matrix(mtx), BCPlaid(), fit.model = ~m + a + b, verbose=T)
Bc
drawHeatmap(as.matrix(mtx),Bc,1, )
drawHeatmap2(as.matrix(mtx), Bc, plotAll=T)
bubbleplot(as.matrix(mtx),Bc,showLabels=T)

modules <- as(Bc, "ISAModules")
# ISA2heatmap(modules, 2, as.matrix(mtx))
# profilePlot(modules, 2, mtx.eset, plot = "samples")
optimalorder <- OrderEV(modules)
optimalorder$status

data<-as.matrix(mtx)
modules <- isa(data)


## Add metadata associated with the rows of the data set
rowdata <- outer(1:nrow(data), 1:sample(1:20, 1), function(x, y) {
  paste("row description (", x, ", ", y, ")", sep="")
})
rownames(rowdata) <- rownames(data)
colnames(rowdata) <- paste("row tag", seq_len(ncol(rowdata)))

## Add metadata associated with the columns of the data set 
coldata <- outer(1:ncol(data), 1:sample(1:20, 1), function(x, y) {
  paste("column description (", x, ", ", y, ")", sep="")
})
rownames(coldata) <- colnames(data)
colnames(coldata) <- paste("column tag", seq_len(ncol(coldata)))

## Merge the different annotations in a single list and 
## add a few global things
description <- list(
  experiment=list(
    title="Title", 
    xaxislabel="x-Axis Label",
    yaxislabel="y-Axis Label",
    name="Author", 
    lab="Address", 
    abstract="Abstract", 
    url="URL", 
    annotation="Annotation", 
    organism="Organism"),
  coldata=coldata,
  rowdata=rowdata
)

ExportEV(modules, data, optimalorder, filename = "F:/file2.evf", description=description)
ExportEV(modules, data, optimalorder, filename="f:/file11.evf", description=description)