# Set up the environment
library(plyr)
library(dplyr)
library(reshape2)
library(cluster)
library(ggplot2)
library(gplots)
library(Hmisc)
library(Biobase)
library(limma)
library(pander)
library(colorRamps)
library(genefilter)
library(xlsx)
source("/Users/mikhail/Documents/Work/GenomeRunner/R.GenomeRunner/genomeRunner_file_formatting_functions.R")
# Work paths
gfAnnot <- read.xlsx2("/Users/mikhail/Documents/Work/GenomeRunner/genomerunner_database/hg19/GFs_hg19_joined_cell_factor.xlsx", sheetName="GFs_hg19_joined_cell_histone_1")
#gfAnnot <- tbl_df(read.table("/Users/mikhail/Documents/Work/GenomeRunner/genomerunner_database/hg19/GFs_hg19_joined_cell_factor.txt", sep="\t", header=F))
#gfAnnot$V1 <- make.names(gfAnnot$V1)
#cellAnnot <- tbl_df(read.table("/Users/mikhail/Documents/Work/GenomeRunner/genomerunner_database/ENCODE_cells.txt", sep="\t", header=T, fill=T, quote="\""))
# Home paths
#gfAnnot <- tbl_df(read.table("/Users/mikhaildozmorov/Documents/Work/GenomeRunner/genomerunner_database/hg19/GFs_hg19_joined_cell_factor.txt", sep="\t", header=T))
# Windows paths
#gfAnnot <- tbl_df(read.table("F:/Work/GenomeRunner/genomerunner_database/hg19/GFs_hg19_joined_cell_factor.txt", sep="\t", header=F))
#cellAnnot <- tbl_df(read.table("F:/Work/GenomeRunner/genomerunner_database/hg19/ENCODE_cells.txt", sep="\t", header=T, fill=T, quote="\""))
cellAnnot <- aggregate(gfAnnot$celldescr, list(gfAnnot$cell), unique)
colnames(cellAnnot) <- c("cell", "description")
cellAnnot <- cellAnnot[ cellAnnot$cell != "", ] # Remove empty cells
# Define color palette
#color<-colorRampPalette(c("blue", "yellow", "red")) # Define color gradient
color <- matlab.like
# Adjust clustering parameters.
# Distance: "euclidean", "maximum","manhattan" or "minkowski". Do not use "canberra" or "binary"
# Clustering: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
granularity <- 10
dist.method <- "euclidean"  
hclust.method <- "ward.D2"
# For coordinate-extracting function
suppressMessages(library(biomaRt))
# mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org",path="/biomart/martservice",archive=FALSE, verbose=TRUE) # Last mart containing HG19 genome annotation
upstream <- 2000 # Definition of the promoter - 2000bp upstream
downstream <- 500 # and 500bp downstream

#' Enrichment analysis data import
#' 
#' A function to load enrichment analysis matrix(es) and remove non-informative epigenomic marks
#'
#' @param dname a string, or a character vector containing one or multiple paths to the matrix file(s). Multiple matrixes will be 'rbind'. If a matrix name has 'PVAL' in its name, the data will be -log10-transformed and filtered to remove rows with nothing significant. If a matrix name has 'OR' in its name, the data will be log2-transformed. If no such keywords are found, the data is returned as is. 
#' @param subset a string used to subset a list of genomic features. Default - none. Examples - "Tfbs", "Histone"
#'
#' @return a matrix of filtered data
#' @export
#' @examples
#' mtx <- load_gr_data("data/ENCODE/matrix_OR.txt")
#' mtx <- load_gr_data(c("data/ENCODE_Tfbs/matrix_PVAL.txt", "data/ENCODE_Histone/matrix_PVAL.txt"), subset=c("Tfbs", "Histone"))
##
load_gr_data <- function(dname, subset="none") {
  # mtx<-do.call("rbind", lapply(dname, function(fn) as.matrix(read.table(paste(fn, "matrix.txt", sep=""), sep="\t", header=T, row.names=1))))
  mtx.list <- list()
  for (d in dname) {
    mtx.list <- c(mtx.list, list(as.matrix(read.table(d, sep="\t", header=T, row.names=1, stringsAsFactors=F, check.names=FALSE))))
  }
  # rbind matching column names
  # https://stackoverflow.com/questions/16962576/how-can-i-rbind-vectors-matching-their-column-names
  mtx <- do.call("rbind", lapply(mtx.list, function(x) x[, match(colnames(mtx.list[[1]]), colnames(x)) ]))
  class(mtx) <- "numeric" # Convert to numbers
  if (subset != "none") {
    mtx <- mtx[grep(paste(subset, collapse="|"), rownames(mtx), ignore.case=T), ] # filter GF list
  }
  # Exploratory: check quantiles and remove diseaases showing no enrichments
  # mtx.sumstat <- as.data.frame(apply(mtx, 2, quantile)) # Get quantiles
  # mtx <- mtx[ , apply(mtx.sumstat, 2, function(x) sum(abs(x)) != 5)] # Remove those that have all "1" or "-1"
  # Optional: filter unused genomic features
  # mtx<-mtx[grep("snp", rownames(mtx), ignore.case=T, invert=T), ]
  if (grepl("PVAL", d)) { # Process matrix of raw p-values
    mtx<-mtx.transform(mtx) # -log10 transform p-values
    # Optional: adjust columns for multiple testing. See utils.R for the function definition.
    # mtx<-mtx.adjust(mtx) 
    dim(mtx) # Check original dimensions
    # Define minimum number of times a row/col should have values above the cutoffs
    numofsig <- 1
    cutoff <- -log10(0.1) # q-value significance cutoff
    # What remains if we remove rows/cols with nothing significant
    dim(mtx[apply(mtx, 1, function(x) sum(abs(x) > cutoff)) >= numofsig, ])
    # apply(mtx, 2, function(x) sum(abs(x)>cutoff))>=numofsig])
    # Trim the matrix
    mtx<-mtx[apply(mtx, 1, function(x) sum(abs(x) > cutoff)) >= numofsig, ]
    # apply(mtx, 2, function(x) sum(abs(x)>cutoff))>=numofsig]
    # If there are columns with SD=0, add jitter to it
    set.seed(1)
    mtx[, apply(mtx, 2, sd) == 0] <- jitter(mtx[, apply(mtx, 2, sd) == 0], factor=0.1)
  }
  if (grepl("OR", d)) { # Prpcess matrix of odds ratios
    mtx <- log2(mtx)
  }
  return(as.matrix(mtx)) # Return (processed) data
}

## ----------------------------------------------------------------------------------
#' Extracts genomic coorfinates for a list of EntrezIDs
#' 
#' A function to extract strand-specific genomic coordinates of a list of EntrezIDs. Homo Sapiens only, for now.
#' Writes the output to fileName
#' 
#' @param entrezIDs a numeric vector of EntrezIDs. Required
#' @param fileName the name of the file to output converted coorcinates
#' 
#' @export
#' 
#' @examples
#' entrez2bed(c("3106","51752","149233","118429","864"), "coords.bed")
##
entrez2bed <- function(entrezIDs, fileName){
  coords <- getBM(attributes=c('chromosome_name','start_position', 'end_position', 'hgnc_symbol', 'strand'), filters='entrezgene', values=entrezIDs, mart=mart, uniqueRows=F)
  coords <- coords[ coords$chromosome_name %in% c(seq(1,22), "X"), ] # Keep only normal chromosomes
  ind.pos <- coords$strand == 1 # Indexes of positive strand
  ind.neg <- coords$strand == -1 # and negative strand
  # Recalculate promoters depending on strand
  coords$end_position[ ind.pos ] <- coords$start_position[ ind.pos ] + downstream
  coords$start_position[ ind.pos ] <- coords$start_position[ ind.pos ] - upstream
  coords$start_position[ ind.neg ] <- coords$end_position[ ind.neg ] - downstream
  coords$end_position[ ind.neg ] <- coords$end_position[ ind.neg ] + upstream
  # Make BED format
  coords$strand[ ind.pos ] <- "+"
  coords$strand[ ind.neg ] <- "-"
  coords$chromosome_name <- paste("chr", coords$chromosome_name, sep="")
  # Format into BED format
  coords <- cbind(coords$chromosome_name, coords$start_position, coords$end_position, coords$hgnc_symbol, rep(".", nrow(coords)), coords$strand)
  # Sort, merge and save the coordinates
  write.table(coords, "tmp.bed", sep="\t", quote=F, row.names=FALSE, col.names=FALSE)
  command <- paste("bedtools sort -i tmp.bed | mergeBed -i - -c 4 -o collapse > ", fileName, sep="")
  system(command)
  unlink("tmp.bed")
}



## ----------------------------------------------------------------------------------
## Convert a matrix of raw p-values (with "-" indicating depletion) into -log10-transformed
mtx.transform<-function(x){
  tmp<- -log10(abs(x)) # -log10 transformation without sign
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (x[i, j]<0) {tmp[i, j]<- -tmp[i,j]} # Add sign, if needed
    }
  }
  return(tmp)
}

## ----------------------------------------------------------------------------------
## Convert a matrix of p-values (with "-" indicating depletion) into Z-scores
mtx.transform.p2z <- function(datapv) {
  dataz = sign(datapv) * qnorm(p = as.matrix(abs(datapv))/2, lower.tail = FALSE)
  return(dataz)
}

## ----------------------------------------------------------------------------------
## Convert a matrix of Z-scores (with "-" indicating depletion) into p-values
mtx.transform.z2p <- function(dataz) {
  datapv = sign(dataz) * 2* pnorm(q = as.matrix(-abs(dataz)))
  return(datapv)
}


## ----------------------------------------------------------------------------------
#' Creates crossproduct correlation matrix
#' 
#' A function to convert matrix of enrichment profiles (columns) into correlation matrix using crossproduct
#' 
#' @param dataz an N x M matrix of signed Z-scores. Required
#' @param quantile top X% of Z-scores to consider for correlations. Default - 5%.
#' 
#' @return an N x N matrix of correlation crossproducts
#' @export
#' 
#' @examples
#' mtx.cr <- mtx.cor.crossprod(dataz)
##
mtx.cor.crossprod <- function(dataz, quantile=0.05) {
  transform01 = function(vec, quantile) {
    threshold = quantile(x = abs(vec), probs = 1-quantile);
    rez = sign(vec) * (abs(vec) >= threshold);
    return( rez );
  }
  data01 = data.frame( lapply(as.data.frame(dataz), transform01, quantile) )
  cr = crossprod(as.matrix(data01));
  return(cr)
}

## ----------------------------------------------------------------------------------
## Convert a matrix of -log10-transformed p-values (with "-" indicating depletion) into raw linear scale
mtx.untransform<-function(x){
  tmp<- 1/10^abs(x) # -log10 transformation without sign
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (x[i, j]<0) {tmp[i, j]<- -tmp[i,j]} # Add sign, if needed
    }
  }
  return(tmp)
}

## ----------------------------------------------------------------------------------
## Correct for multiple testing the matrix of transformed p-values
mtx.adjust.log10<-function(x, method="fdr"){ 
  tmp<- -log10(apply(1/10^abs(x), 2, p.adjust, method=method)) # Adjust absolute p-values for multiple testing and convert them back to -log10
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (x[i,j]<0) {tmp[i,j]<- -1*tmp[i,j]} # Add sign, if needed
    }
  }
  return(tmp)
}

## ----------------------------------------------------------------------------------
## Adjust raw p-values for mutliple testing
mtx.adjust.raw <- function(x, adjust="fdr") {
  tmp <- x; # Keep signs in the original version
  tmp.adj <- apply(abs(x), 2, p.adjust, method=adjust) 
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (tmp[i, j] < 0) tmp.adj[i, j] <- -1*tmp.adj[i, j] # Add sign, if needed
    }
  }
  return(tmp.adj)
}

## ----------------------------------------------------------------------------------
## Adjust ONE vector of p-values for mutliple testing, return -log10-transformed values
mtx.adjust.1 <- function(x, adjust="fdr", isLog10=TRUE) { 
  if (isLog10) {
    tmp <- -log10(p.adjust(1/10^abs(x), method=adjust))
  } else {
    tmp <- -log10(p.adjust(abs(x), method=adjust))
  }
  for (i in 1:length(x)) {
    if (x[i] < 0) tmp[i] <- -1*tmp[i] # Add a "-" sign, if it was there previously  
  }
  return(tmp)
}

## ----------------------------------------------------------------------------------
## Filters non-significant rows from a matrix
mtx.filter <- function(x, pval=0.05){
  idx <- apply(x, 1, function(i) {sum(abs(i) < pval) >= 1})
  tmp <- as.matrix(x[idx, , drop=F])
  colnames(tmp) <- colnames(x)
  #rownames(tmp) <- rownames(x)[idx]
  idxEnrich <- as.logical(apply(tmp, 1, function(x) {length(x[x > 0 & x < pval]) > 0}))
  tmp1 <- as.matrix(tmp[idxEnrich, , drop=F])
  colnames(tmp1) <- colnames(tmp)
  #rownames(tmp1) <- rownames(tmp)[idxEnrich]
  idxDeplet <- as.logical(apply(tmp, 1, function(x) {length(x[x < 0 & -x < pval]) > 0}))
  tmp2 <- as.matrix(tmp[idxDeplet, , drop=F])
  colnames(tmp2) <- colnames(tmp)
  #rownames(tmp2) <- rownames(tmp)[idxDeplet]
  if (nrow(tmp1) > 1) tmp1 <- tmp1[order(apply(tmp1, 1, function(x) {mean(x[x > 0 & x < pval])}), decreasing=F), , drop=F]
  if (nrow(tmp2) > 1) tmp2 <- tmp2[order(apply(tmp2, 1, function(x) {mean(x[x < 0 & -x < pval])}), decreasing=T), , drop=F]
  return(list(tmp1 <- data.frame(gfs=rownames(tmp1), tmp1), 
              tmp2 <- data.frame(gfs=rownames(tmp2), tmp2)))
}

## ----------------------------------------------------------------------------------
## Shows top 10 enriched and depleted associations as tables
##
## Parameters:
##     fname - full path to enrichment matrix
##     isLog10 - whether a matrix is log10-transformed (TRUE for GR, FALSE for GR WEB)
##     adjust - whether to correct p-values for multiple testing ("fdr" for GR, "none" for GR WEB)
##     pval - cutoff of calling an enrichment significant
##
showTable <- function(fname, isLog10=TRUE, adjust="fdr", pval=0.05) {
  mtx <- read.table(fname, sep="\t", row.names=1, header=T)
  if (isLog10) {mtx <- mtx.untransform(mtx)} # Un- -log10 transform p-values, if needed
  mtx <- mtx.adjust.raw(mtx, adjust) # Adjust for multiple testing
  mtx <- mtx.filter(mtx, pval) # Get two sepatate lists, enchched and depleted
  numEnrich <- nrow(mtx[[1]]) # Total number of significantly enriched
  numDeplet <- nrow(mtx[[2]]) # and depleted associations
  if (numEnrich > 0) {
    print(paste("The total number of significantly ENRICHED associations is:", numEnrich))
    if (numEnrich > 10) numEnrich <- 10
    #mtx.enrich <- merge(data.frame(pval=format(as.matrix(mtx[[1]][1:numEnrich, , drop=F]), digits=3, scientific=T),
    #                               row.names=rownames(as.matrix(mtx[[1]]))[1:numEnrich]), trackDb.hg19, by="row.names", all.x=T, sort=T)
    #mtx.enrich <- data.frame(pval=format(as.matrix(mtx[[1]][1:numEnrich, , drop=F]), digits=3, scientific=T), row.names=rownames(as.matrix(mtx[[1]]))[1:numEnrich])
    mtx.enrich <- left_join(mtx[[1]], trackDb.hg19, by=c("gfs" = "V1"))
    pander(mtx.enrich[1:numEnrich, ], split.table=Inf)
    mtx[[1]] <- as.data.frame(mtx[[1]]); rownames(mtx[[1]]) <- mtx[[1]][, 1]; mtx[[1]] <- mtx[[1]] <- mtx[[1]][, -1, drop=F]
    mtx[[1]][mtx[[1]] == 0] <- .Machine$double.xmin
    barplot1(mtx.transform(mtx[[1]][1:numEnrich, , drop=F]), 15)
    #    grid.table(mtx.enrich, gp=gpar(fontsize=6))
    print("---------------------------------------------------------------")
  }
  if (numDeplet > 0) {
    print(paste("The total number of significantly DEPLETED associations is:", numDeplet))
    if (numDeplet > 10) numDeplet <- 10
    mtx.deplet <- merge(data.frame(pval=format(as.matrix(mtx[[2]][1:numDeplet, , drop=F]), digits=3, scientific=T), row.names=rownames(as.matrix(mtx[[2]]))[1:numDeplet]), trackDb.hg19, by="row.names", all.x=T)
    #mtx.deplet <- data.frame(pval=format(as.matrix(mtx[[2]][1:numDeplet, , drop=F]), digits=3, scientific=T), row.names=rownames(as.matrix(mtx[[2]]))[1:numDeplet])
    mtx.deplet <- left_join(mtx[[2]], trackDb.hg19, by=c("gfs" = "V1"))
    pander(mtx.deplet[1:numDeplet, ], split.table=Inf)
    mtx[[2]] <- as.data.frame(mtx[[2]]); rownames(mtx[[2]]) <- mtx[[2]][, 1]; mtx[[2]] <- mtx[[2]] <- mtx[[2]][, -1, drop=F]
    mtx[[2]][mtx[[2]] == 0] <- .Machine$double.xmin
    barplot1(mtx.transform(mtx[[2]][1:numDeplet, , drop=F]), 15)
    #grid.table(mtx)
    print("---------------------------------------------------------------")
  }
}

## ----------------------------------------------------------------------------------
## Make a barplot of -log10-transformed p-values matrix.
## 
## Parameters:
##     mtx - matrix of -log10-transformed p-values
##     location - where to put the legend (e.g., "topright", "bottomright")
##     bottom - bottom margin (e.g., 5 or 15)
##     names.args - names of bar groups
##     pval - where to draw cutoff lines
## 
barplot1<-function(mtx, location="topright", bottom=5, names.args, pval=0.1){
  par(mar=c(bottom,5,2,2)+0.1)
  groupcolors<-rainbow(ncol(mtx)) #c("yellow2","steelblue3","steelblue3","springgreen)
  if (grepl("top", location)) {
    txt <- "Overrepresented regulatory associations"
  } else {
    txt <- "Underrepresented regulatory associations"
  }
  mtx[mtx == Inf] <- 308 # Replace infinite values to a finite number
  b<-barplot(as.matrix(t(mtx)), beside=T,  ylab="-log10(p-value)\nnegative = underrepresentation", col=groupcolors,space=c(0.2,1), cex.names=0.6, las=2, names.arg=names.args, main=txt) # ,legend.text=colnames(mtx),args.legend=list(x=7,y=4))
  lines(c(0,100),c(-log10(pval),-log10(pval)),type="l",lty="dashed",lwd=2)
  lines(c(0,100),c(log10(pval),log10(pval)),type="l",lty="dashed",lwd=2)
  legend(location, legend=colnames(mtx), fill=groupcolors, cex=0.6)
 }

## ----------------------------------------------------------------------------------
## Trim the rows/columns of the enrichment matrix until each row/column has at least
## numofsig number of cells above pval threshold
##
mtx.trim.numofsig <- function(mtx.cast, pval=0.1, numofsig=1) {
dims.new <- dim(mtx.cast[apply(mtx.cast, 1, function(x) sum(abs(x) > -log10(pval), na.rm=T)) >= numofsig, 
                         apply(mtx.cast, 2, function(x) sum(abs(x) > -log10(pval), na.rm=T)) >= numofsig])
repeat {
    # Trim the matrix
    mtx.cast<-mtx.cast[apply(mtx.cast, 1, function(x) sum(abs(x) > -log10(pval), na.rm=T)) >= numofsig, 
                       apply(mtx.cast, 2, function(x) sum(abs(x) > -log10(pval), na.rm=T)) >= numofsig] 
    dims.old <- dim(mtx.cast)
    dims.new <- dim(mtx.cast[apply(mtx.cast, 1, function(x) sum(abs(x) > -log10(pval), na.rm=T)) >= numofsig, 
                             apply(mtx.cast, 2, function(x) sum(abs(x) > -log10(pval), na.rm=T)) >= numofsig])
    if (all(dims.new == dims.old)) {
      return(mtx.cast)
      break
    }
  }
}

## ----------------------------------------------------------------------------------
## Trim the rows/columns of the enrichment matrix until each row/column has less NAs than
## the numofnas number
##
mtx.trim.numofnas <- function(mtx.cast, numofnas=1) {
  dims.new <- dim(mtx.cast[apply(mtx.cast, 1, function(x) sum(is.na(x))) <= numofnas,
                           apply(mtx.cast, 2, function(x) sum(is.na(x))) <= numofnas])
  repeat {
    # Trim the matrix
    mtx.cast<-mtx.cast[apply(mtx.cast, 1, function(x) sum(is.na(x))) <= numofnas,
                       apply(mtx.cast, 2, function(x) sum(is.na(x))) <= numofnas]
    dims.old <- dim(mtx.cast)
    dims.new <- dim(mtx.cast[apply(mtx.cast, 1, function(x) sum(is.na(x))) <= numofnas, 
                             apply(mtx.cast, 2, function(x) sum(is.na(x))) <= numofnas])
    if (all(dims.new == dims.old)) {
      return(mtx.cast)
      break
    }
  }
}

## ----------------------------------------------------------------------------------
## Creates Cell x Factor heatmap from a Histone/TFBS one-column matrix of enrichments.
## Also, plots barplots of most enriched/depleted associations
##
## Parameters:
##     fname - full path to enrichment matrix
##     colnum - which column(s) to use, in case of multi-column matrix. If a single column
## selected, and the factor is set, a heatmap of Cell x Factor enrichments is plotted
##     factor - subset matrix by "Histone"/"Tfbs" enrichments, or "none". If "none", 
## barplot is labeled with table names, otherwise, by "cell.factor" labels. Multiple factors are 'AND'
## allowed, e.g., using c("Histone", "Gm12878") will select histone AND Gm12878 datasets
##     cell - subset matrix by cell type(s). Multiple terms are OR allowed, e.g., using 
## c("Gm12878", "K562") will select Gm12878 OR K562 datasets
##     isLog10 - whether a matrix is log10-transformed (TRUE for GR, FALSE for GR WEB). A special case is when
## isLog10="OR" and a matrix of odds ratios is provided - it will be log2-transformed. For odds ratios, no adjustment
## for multiple testing is done, and pval/numtofilt parameters should be 1 and 1, respectively. Mind what's going on
##     adjust - whether to correct p-values for multiple testing ("fdr" for GR, "none" for GR WEB)
##     pval - cutoff of calling an enrichment significant
##     numtofilt - number of cells in each row/column to filter. Valid for 'heat' plot only
##     toPlot - What to plot. Can be "heat", "bar" or "barup"/"bardn", "lines", "corrPearson" or "corrSpearman". 
## "heat" is relevant for 1-column Histone/Tfbs/BroadHmm-specific matrix. Returns clustered carpet
## "bar" can be used on multiple columns, matrix can be filtered by factor/cell. Returns either up- or downregulated matrixes, or a list of both
## "lines" used to compare profiles of p-values across multiple columns, matrix can be filtered. Returns nothing
## "corrPearson"/"corrSpearman" - correlation heatmap. Works best on multi-column/row matrix. Returns clustered carpet
##     fileName - save the reshaped wide matrix into a file
##
## Examples:
##     For the original GR:
## showHeatmap("matrix.txt", colnum=1, factor="Histone", cell="none", isLog10=TRUE, adjust="fdr", pval=0.1, numtofilt=4, toPlot="heat", fileName=NULL)
##     For the GR WEB:
## showHeatmap("matrix.txt", colnum=1, factor="Tfbs", cell="none", isLog10=FALSE, adjust="none", pval=0.1, numtofilt=4, toPlot="heat", fileName=NULL)
## showHeatmap("matrix.txt", colnum=c(1, 5), factor="none", cell="none", isLog10=FALSE, adjust="none", pval=0.1, numtofilt=1, toPlot="bar", fileName=NULL)
## showHeatmap("matrix.txt", colnum=c(5,6,7,8), factor=c("Tfbs"), cell="none", isLog10=FALSE, adjust="none", pval=0.1, numtofilt=6, toPlot="lines")
## showHeatmap("matrix.txt", colnum=seq(1,50), factor="none", cell="none", isLog10=FALSE, adjust="none", pval=0.5, numtofilt=1, toPlot="corrSpearman")
##
## colnum=1; factor="none"; cell="none"; isLog10=FALSE; adjust="fdr"; pval=0.1; numtofilt=1; toPlot="bar"; fileName=NULL

showHeatmap <- function(fname, colnum=1, factor="none", cell="none", isLog10=FALSE, adjust="fdr", pval=0.1, numtofilt=1, toPlot="bar", fileName=NULL) {
  mtx <- tbl_df(read.table(fname, sep="\t", fill=T, header=F, stringsAsFactors=F)) # No header and row.names
  cols <- mtx[1, ] # Keep header
  # Subsetting by factor
  if (factor != "none") {
    mtx <- mtx[grep(factor[1], mtx$V1), c(1, colnum + 1)] # Subset by factor and column 
    if (length(factor) > 1) {
      for (i in 2:length(factor)) {
        mtx <- mtx[grep(factor[i], mtx$V1), ] # Subset by factor, keep columns
      }  
    }
  } else {
    mtx <- mtx[-1, c(1, colnum + 1)] # Do not subset
  }
  # Subsetting by cell
  if (cell != "none") {
    mtx <- mtx[grepl(paste(cell, collapse="|"), mtx$V1), ]
  }
  # Adjust for multiple testing and -log10 transform, if needed
  for (i in 1:length(colnum)) {
    if(is.logical(isLog10)){ # If p-values are provided, which is defined by TRUE or FALSE isLog10, adjust accordingly
      mtx[, i + 1] <- mtx.adjust.1(as.numeric(unlist(mtx[, i + 1])), adjust=adjust, isLog10=isLog10)
    } else { # If odds ratios, which is defined by isLog10 equal to non-logical value like "OR", simply keep the numerical values. Remember OR ranges 0-1 for underrepresentation and 1-Inf for overrepresentation
      mtx[, i + 1] <- as.numeric(unlist(mtx[, i + 1]))
    }
  }
  # Join with annotations
  mtx <- left_join(mtx, gfAnnot[, c("name", "cell", "factor", "description")], by=c("V1" = "name"))
  # Assign columns
  colnames(mtx) <- c("GF", make.names(cols[colnum ]), "cell", "factor", "description") # Rename columns  
  # Save the matrix, if needed
  if (!is.null(fileName)) { 
    write.table(mtx, fileName, sep="\t", row.names=F)
  }
  
  ## Creates Cell x Factor heatmap from a matrix of enrichments from a Histone/Tfbs matrix.
  if (length(colnum) == 1 & (toPlot == "heat")) { # If only 1 column selected, we can plot heatmap
    # Make wide matrix. 
    pmax <- function(x) { x[order(abs(x), decreasing=T)][1] } # Get absolute maximum p-value, keeping sign
    mtx.cast <- dcast(mtx, formula=cell~factor, fun.aggregate=pmax, value.var=make.names(cols[colnum]))
    
    # Make wide matrix. To properly handle duplicates, use https://stackoverflow.com/questions/12831524/can-dcast-be-used-without-an-aggregate-function
#     tmp1 <- ddply(mtx, .(cell, factor), transform, newid = paste(cell, seq_along(factor)))
#     out <- dcast(tmp1, cell + newid ~ factor, value.var=make.names(cols[colnum]))
#     out <- out[,-which(colnames(out) == "newid")]
#     mtx.cast <- out; rm(tmp1, out)

    rownames(mtx.cast) <- make.names(mtx.cast$cell, unique=T) # Reassign rownames
    mtx.cast <- mtx.cast[, -1] # Remove no longer needed first column
    mtx.cast <- mtx.trim.numofsig(mtx.cast, pval=pval, numofsig=numtofilt) # Filter by counting number of significant cells
    # mtx.cast <- mtx.trim.numofnas(mtx.cast, numofnas=numtofilt) # Not working currently
    if (nrow(mtx.cast) == 0 | ncol(mtx.cast) == 0) { 
      print("Nothing significant, cannot plot heatmap")
      return()
    } else {
      # Plotting
      par(cex.main=0.65, oma=c(5,0,0,5), mar=c(5, 4.1, 4.1, 5)) # Adjust margins
      dist.method<-"euclidean"  
      hclust.method<-"ward.D2"
      notecex <- -0.05*ncol(mtx.cast) + 1 # Size of cell text
      if (notecex < 0.5) { notecex <- 0.5 }
      mtx.plot <- as.matrix(mtx.cast) # Matrix to plot
      mtx.max <- max(abs(mtx.plot[mtx.plot != max(abs(mtx.plot), na.rm=T)]), na.rm=T) # Second to max value
      my.breaks <- c(seq(-mtx.max, 0, length.out=10), 0, seq(0, mtx.max, length.out=10)) # Breaks symmetric around 0
      h<-heatmap.2(mtx.plot, trace="none", density.info="none", col=color, distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, 
                   cexRow=0.8, cexCol=0.8, breaks=my.breaks,
                   cellnote=formatC(1/10^abs(as.matrix(mtx.cast)), format="e", digits=2), notecol="black", notecex=notecex)
      return(h$carpet)
    }    
  } else if (grepl("bar", toPlot)) { # If more than 1 column, we can also plot barplots
    ## Plot barplot representation of the enrichments
    mtx <- as.data.frame(mtx) # Make data frame, to allow row names
    rownames(mtx) <- mtx$GF; mtx <- mtx[, -1] # Make row names
    # Summarize the values by cell type and factor
    pmax <- function(x) { x[order(abs(x), decreasing=T)][1] } # Get absolute maximum p-value, keeping sign
    mtx <- mtx %>% dplyr::select(-description) %>% group_by(cell, factor) 
    mtx <- melt(mtx, id.vars=c("cell", "factor"), , variable.name = "experiment", value.name = "pvalue")
    mtx <- mtx %>% group_by(cell, factor, experiment) %>% dplyr::summarise(pmax = pmax(pvalue))
    mtx <- dcast(mtx, cell + factor ~ experiment, value.var="pmax")
    mtx <- cbind(dplyr::select(mtx, 3:ncol(mtx)), dplyr::select(mtx, cell, factor))
    # Get sorted lists
    mtx.sorted.up <- list(); mtx.sorted.dn <- list() # Storage for sorted matrixes 
    for (i in 1:length(colnum)) {
      # For each column, reorder p-values and store top X most significantly enriched
      mtx.sorted.up[[length(mtx.sorted.up) + 1]] <- mtx[order(mtx[, i], decreasing=T)[1:round(30/length(colnum))], ]
      # And depleted
      mtx.sorted.dn[[length(mtx.sorted.dn) + 1]] <- mtx[order(mtx[, i], decreasing=F)[1:round(30/length(colnum))], ]
    }
    # Combine lists into matrixes
    mtx.barplot.up <- ldply(mtx.sorted.up, rbind)
    mtx.barplot.dn <- ldply(mtx.sorted.dn, rbind)
    # Remove values going in wrong directions
    mtx.barplot.up <- mtx.barplot.up[ apply(mtx.barplot.up[, seq(1:length(colnum)), drop=F], 1, function(x) { any(x > -log10(pval)) }), , drop=F]
    mtx.barplot.dn <- mtx.barplot.dn[ apply(mtx.barplot.dn[, seq(1:length(colnum)), drop=F], 1, function(x) { any(x < log10(pval)) }), , drop=F]
    # Remove rows with NAs
    mtx.barplot.up <- mtx.barplot.up[ complete.cases(mtx.barplot.up), , drop=F]
    mtx.barplot.dn <- mtx.barplot.dn[ complete.cases(mtx.barplot.dn), , drop=F]
    # Reassign names
    names.args.up <- paste(mtx.barplot.up$cell, mtx.barplot.up$factor, sep=":")
    names.args.dn <- paste(mtx.barplot.dn$cell, mtx.barplot.dn$factor, sep=":")
    #names.args.up[names.args.up == "NA:NA"] <- make.names(unlist(lapply(mtx.sorted.up, rownames)), unique=T)[names.args.up == "NA:NA"]
    #names.args.dn[names.args.dn == "NA:NA"] <- make.names(unlist(lapply(mtx.sorted.dn, rownames)), unique=T)[names.args.dn == "NA:NA"]
    bottom <- 8
    # Plot barplots
    if (!grepl("dn", toPlot) & (nrow(mtx.barplot.up) > 0)) { barplot1(mtx.barplot.up[, seq(1:length(colnum)), drop=F], "topright", bottom=bottom, names.args=names.args.up, pval=pval) }
    if (!grepl("up", toPlot) & (nrow(mtx.barplot.dn) > 0 )) { barplot1(mtx.barplot.dn[, seq(1:length(colnum)), drop=F], "bottomright", bottom=bottom, names.args=names.args.dn, pval=pval) }
    if ((!grepl("up", toPlot)) & (!grepl("dn", toPlot))) {
      return(list(up=mtx.barplot.up, dn=mtx.barplot.dn))
    } else if (grepl("up", toPlot)) {
      return(mtx.barplot.up)
    } else {
      return(mtx.barplot.dn)
    }
  } else if ((length(colnum) > 1) & (toPlot == "lines")) {
    ## Plot lines. http://kohske.wordpress.com/2010/12/27/faq-geom_line-doesnt-draw-lines/
    df <- mtx[, 1:(length(colnum) + 1)]
    df$GF <- factor(df$GF)
    df.melted <- melt(df, id.vars="GF")
    ggplot(df.melted, aes(x=variable, y=value, colour=GF, group=GF)) + geom_line() + geom_point() + theme(legend.position="none")
  } else if ((length(colnum) > 3) & (grepl("corr", toPlot))) {
    ## Plot correlation heatmap from the original multi-column matrix
    mtx <- as.data.frame(mtx) # Make data frame, to allow row names
    rownames(mtx) <- mtx$GF; mtx <- mtx[, -1] # Make row names
    mtx <- mtx.trim.numofsig(mtx[, 1:length(colnum) ], pval=pval, numofsig=numtofilt) # Filter by counting number of significant cells
    if (grepl("Pearson", toPlot)) {
      mtx.cor <- rcorr(as.matrix(mtx), type="pearson")
    } else {
      mtx.cor <- rcorr(as.matrix(mtx), type="spearman")
    }
    par(cex.main=0.65, oma=c(5,0,0,5), mar=c(5, 4.1, 4.1, 5)) # Adjust margins
    color<-colorRampPalette(c("blue", "yellow")) # Define color gradient
    dist.method<-"euclidean"  
    hclust.method<-"ward.D2"
    notecex <- -0.05*ncol(mtx.cor[[1]]) + 1 # Size of cell text
    if (notecex < 0.4) { notecex <- 0.4 }
    granularity <- 7
    my.breaks <- seq(min(mtx.cor[[1]][mtx.cor[[1]]!=min(mtx.cor[[1]])]),
                     max(mtx.cor[[1]][mtx.cor[[1]]!=max(mtx.cor[[1]])]),
                     length.out=(2*granularity + 1))
    h<-heatmap.2(as.matrix(mtx.cor[[1]]), trace="none", density.info="none", symkey=T, col=color, distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, 
                 cexRow=0.8, cexCol=0.8, breaks=my.breaks, cellnote=formatC(as.matrix(mtx.cor[[1]]), format="f", digits=2), notecol="black", notecex=notecex)
    return(h$carpet)
    }
}

## ----------------------------------------------------------------------------------
## Summarizes a matrix of n GFs x m FOIs by cell-factor most significant p-value
## The rationale here is that different institutions provide data for the same cell and factor,
## e.g., H1hesc:H3K4me3 from UW and H1hesc:H3K4me3 from Broad. The function selects the most
## significant p-value for H1hesc:H3K4me3, and removes duplicated measurements
## Also, returns the summarized matrix
##
## Parameters:
##     mtx - n X m matrix. GF names (rows, n) should be able to be joined with a GF annotation table
##     factor - subset matrix by "Histone"/"Tfbs" enrichments, or "none". Multiple factors 
## are 'OR' allowed, e.g., using c("Histone", "Gm12878") will select histone OR Gm12878 GFs
##     cell - subset matrix by cell type(s). Multiple terms are OR allowed, e.g., using 
## c("Gm12878", "K562") will select Gm12878 OR K562 datasets
##     fileName - save the reshaped wide matrix into a file
##
## Example:
## mtx.summarize(mtx, factor="Histone", cell=c("H1hesc", "H7es", "H9es"), fName="results/summary_all_embryo_hist.txt")
##

mtx.summarize <- function(mtx, factor="none", cell="none", fName=NULL) {
  mtx.annot <- left_join(data.frame(V1=rownames(mtx), mtx), gfAnnot[, c("name", "cell", "factor")], by=c("name" = "V1")) 
  if (factor != "none") { 
    mtx.annot <- mtx.annot %>% dplyr::filter(grepl(paste(factor, collapse="|"), V1)) 
  }
  if (cell != "none") {
    mtx.annot <- mtx.annot %>% dplyr::filter(grepl(paste(cells, collapse="|"), V1))
  }
  mtx.annot <- mtx.annot %>%  dplyr::select(-V1)
  mtx.annot <- melt(mtx.annot, id.vars = c("V3", "V5"), variable.name = "experiment", value.name = "pvalue")
  pmax <- function(x) { x[order(abs(x), decreasing=T)][1] } # Get absolute maximum p-value, keeping sign
  mtx.annot <- mtx.annot %>% group_by(V5, experiment) %>% dplyr::summarise(mean = pmax(pvalue))
  colnames(mtx.annot) <- c("factor", "experiment", "pmax")
  mtx.annot <- dcast(mtx.annot, experiment ~ factor, value.var="pmax")
  rownames(mtx.annot) <- mtx.annot[, 1]
  mtx.annot <- mtx.annot[, -1]
  mtx.annot <- mtx.untransform(as.matrix(mtx.annot))
  if (!is.null(fName)) {
    write.table(mtx.annot, fName, sep="\t", col.names=NA, quote=F)
  }
  return(mtx.annot)
}
