source("utils.R")
library(Hmisc) # For rcorr function
library(gplots)
library(Biobase)
library(limma)

# Read in matrix.txt
dname<-"data//more15_vs_tfbsEncode//" # Which data subfolder to use, change to analyze different data
mtx<-as.matrix(read.table(paste(dname, "matrix.txt", sep=""), sep="\t", header=T, row.names=1))
# dim(na.omit(mtx)) # Just in case, check if NAs are present
# mtx<-na.omit(mtx) # Remove NA rows
mtx<-mtx.transform(mtx) # -log10 transform p-values
# mtx<-mtx.adjust(mtx) # Optional: adjust for multiple testing

# Remove non-significant rows/columns from the matrix
dim(mtx) # Check original dimensions
cutoff<-2 # raw p-value significance cutoff is 0.01. Use 1.3 to have 0.05 significance cutoff
dim(mtx[apply(mtx, 1, function(x){sum(abs(x)>cutoff)})>0, apply(mtx, 2, function(x){sum(abs(x)>cutoff)})>0]) # What remails if remove rows/cols with nothing significant
mtx<-mtx[apply(mtx, 1, function(x){sum(abs(x)>cutoff)})>0, apply(mtx, 2, function(x){sum(abs(x)>cutoff)})>0] # Do it