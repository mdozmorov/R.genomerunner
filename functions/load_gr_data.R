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
  # Load matrix(es) from a (vector of) file(s) located at dname
  mtx.list <- list()
  for (d in dname) {
    mtx.list <- c(mtx.list, list(as.matrix(read.table(d, sep="\t", header=T, row.names=1, stringsAsFactors=F, check.names=FALSE), drop=FALSE)))
  }
  # rbind matrixes, matching column names
  # https://stackoverflow.com/questions/16962576/how-can-i-rbind-vectors-matching-their-column-names
  mtx <- do.call("rbind", lapply(mtx.list, function(x) x[, match(colnames(mtx.list[[1]]), colnames(x)), drop = FALSE ]))
  class(mtx) <- "numeric" # Convert to numbers
  # filter GF list, if specified
  if (subset != "none") {
    mtx <- mtx[grep(paste(subset, collapse="|"), rownames(mtx), ignore.case=T), ] 
  }
  # Transform PVAL and OR matrixes accordingly
  if (grepl("PVAL", d)) { 
    mtx <- mtx.transform(mtx) # -log10 transform p-values
  }
  if (grepl("OR", d)) {
    mtx <- log2(mtx) # log2 transform odds ratios
  }
  # Trim the matrix
  mtx <- mtx[ apply(mtx, 1, function(x) sum(!is.na(x))) > 0, apply(mtx, 2, function(x) sum(!is.na(x))) > 0, drop=FALSE] # Remove rows/columns with all NAs
  mtx <- mtx[ !(apply(mtx, 1, function(x) sum(x == 0) == ncol(mtx))), , drop=F] # If all values in a row are 0, remove these rows
  
  return(as.matrix(mtx)) # Return (processed) data
}