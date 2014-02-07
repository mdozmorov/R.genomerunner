## ----loadLibraries, echo=TRUE, cache=FALSE-------------------------------
source("utils.R")
suppressMessages(library(Hmisc)) # For rcorr function
suppressMessages(library(gplots))
suppressMessages(library(Biobase))
suppressMessages(library(limma))


## ----loadData, echo=c(-5, -6)--------------------------------------------
# Define output and data subfolders to use, change to analyze different data
rname<-"data//results//" # Output folder
# One or more GenomeRunner Web results data folders
dname<-c("data//more15_vs_tfbsEncode//", "data//more15_vs_Other//", "data//more15_vs_ENCODE//")
mtx<-do.call("rbind", lapply(dname, function(fn) as.matrix(read.table(paste(fn, "matrix.txt", sep=""), sep="\t", header=T, row.names=1))))
# Optional: filter unused genomic features
# mtx<-mtx[grep("snp", rownames(mtx), ignore.case=T, invert=T), ]
mtx<-mtx.transform(mtx) # -log10 transform p-values
# Optional: adjust columns for multiple testing. See utils.R for the function definition.
mtx<-mtx.adjust(mtx) 


## ----preprocessData, echo=TRUE, cache=TRUE, dependson='loadData'---------
dim(mtx) # Check original dimensions
# Define minimum number of times a row/col should have values above the cutoffs
numofsig<-1
cutoff<- -log10(0.01) # p-value significance cutoff is 0.01
# What remains if we remove rows/cols with nothing significant
dim(mtx[apply(mtx, 1, function(x) sum(abs(x)>cutoff))>=numofsig, 
        apply(mtx, 2, function(x) sum(abs(x)>cutoff))>=numofsig])
# Trim the matrix
mtx<-mtx[apply(mtx, 1, function(x) sum(abs(x)>cutoff))>=numofsig, 
         apply(mtx, 2, function(x) sum(abs(x)>cutoff))>=numofsig]