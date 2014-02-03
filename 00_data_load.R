## ----loadLibraries, echo=TRUE, cache=FALSE-------------------------------
source("utils.R")
suppressMessages(library(Hmisc)) # For rcorr function
suppressMessages(library(gplots))
suppressMessages(library(Biobase))
suppressMessages(library(limma))


## ----loadData, echo=c(-5, -6)--------------------------------------------
# Define data subfolder to use, change to analyze different data
dname<-"data//more15_vs_tfbsEncode//"
mtx<-as.matrix(read.table(paste(dname, "matrix.txt", sep=""), sep="\t", header=T, row.names=1))
mtx<-mtx.transform(mtx) # -log10 transform p-values
# Optional: adjust columns for multiple testing. See utils.R for the function definition.
# mtx<-mtx.adjust(mtx) 


## ----preprocessData, echo=TRUE, cache=TRUE, dependson='loadData'---------
dim(mtx) # Check original dimensions
cutoff<-2 # raw p-value significance cutoff is 0.01
# What remains if we remove rows/cols with nothing significant
dim(mtx[apply(mtx, 1, function(x){sum(abs(x)>cutoff)})>0, 
        apply(mtx, 2, function(x){sum(abs(x)>cutoff)})>0])
# Trim the matrix
mtx<-mtx[apply(mtx, 1, function(x){sum(abs(x)>cutoff)})>0, 
         apply(mtx, 2, function(x){sum(abs(x)>cutoff)})>0]