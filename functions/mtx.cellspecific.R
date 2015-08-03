#' Cell type-specific enrichment analysis
#' 
#' A function to perform cell type-specific enrichment analysis on a matrix of the enrichment results (raw p-values). This version compares distributions of the overall and cell-specific -log10 p-values using Wilcoxon test.
#' 
#' @param mtx a matrix of the enrichment results (transformed p-values). The cell type-specific enrichment analysis is performed on each column (for each SNP set).
#' @param fname a path to a filename to save the results. Should be with "xlsx" extension. The results for each column are saved in separate worksheets. Default - NULL, do not save.
#' @param pval a p-value cutoff for cell type enrichment significance. Default: 0.01
#' @return A list of the results
#' @export
#' @examples
#' mtx.cellspecific(mtx)
##
mtx.cellspecific <- function(mtx, fname=NULL, pval=0.01) {
  if (nrow(mtx) <= 5) return("Insufficient data for performing cell type-specific enrichment analysis") # If too few genomic features, no analysis can be performed
  n.fois <- ncol(mtx) # Total number of SNP sets (FOIs) to calculate the cell type-specific p-values
  # Prepare the matrix for merging with GF annotations
  mtx <- as.data.frame(cbind(GF=rownames(mtx), mtx))
  mtx <- left_join(mtx, gfAnnot[, c("file_name", "cell", "factor")], by=c("GF" = "file_name")) 
  mtx$cell[ is.na(mtx$cell) ] <- "dummy_cell" # If some file names is not in the gfAnnot dataframe (e.g., user-provided data), 'cell' column will contain NAs. replace them with dummy text to allow cell type-specific analysis
  
  cells <- unique(mtx$cell) # All cell types
  # Global counts
  tot.tests <- nrow(mtx) # Total number of enrichment analyses
  cells.tests <- vector(mode="numeric", length=length(cells)) # Number of analyses per cell type
  names(cells.tests) <- cells
  for(c in 1:length(cells)) {
    cells.tests[c] <- length(mtx$cell[ mtx$cell == cells[c]]) # Number of analyses per cell type
  }
  # Keep cell types that have at least 5 measures
  cells.tests <- cells.tests[ cells.tests > 5 ]
  cells <- names(cells.tests) # Cell types that have sufficient data
  if (length(cells.tests) <= 1) return("Insufficient data for performing cell type-specific enrichment analysis") # If no cells have >5 measures, or if there's only one cell type, no analyses can be done
  
  
  # Column-specific counts
  pval.foi <- list() # for p-values
  stats.foi <- list() # for 2x2 tables
  for(i in 2:(n.fois + 1)){ # Columns are now shifted by 1
    tot.sig <- as.numeric(mtx[, i])  # Column-specific total number of significant results
    cells.sig <- vector(mode="list", length=length(cells)) # Column-specific & cell type-specific number of significant results
    for(c in 1:length(cells)) {
      mtx.sel <- as.numeric(mtx[ mtx$cell == cells[c], i ]) # Column- and cell type-specific vector
      cells.sig[c] <- list(mtx.sel) # How many are significant
    }
    
    pval.foi.cell <- vector(mode="numeric", length=length(cells)) # A vector to store disease- and cell type-specific p-values
    stats.foi.cell <- list() # A list to store disease- and cell type-specific 2x2 tables
    for(c in 1:length(cells)) {
        cells.test <- wilcox.test(cells.sig[[c]], tot.sig, alternative = "greater")
        pval.foi.cell[c] <- cells.test$p.value # Store the enrichment p-values
        stats.foi.cell[c] <- list(c(cells.tests[c], mean(cells.sig[[c]]), mean(tot.sig)))
    }
    names(pval.foi.cell) <- cells # Name the collected vectors
    names(stats.foi.cell) <- cells # as cell names
    pval.foi <- c(pval.foi, list(pval.foi.cell)) # Store them
    stats.foi <- c(stats.foi, list(stats.foi.cell)) # in disease-specific lists
  }
  names(pval.foi) <- colnames(mtx)[2:(n.fois + 1)] # Name the disease-specific lists
  names(stats.foi) <- colnames(mtx)[2:(n.fois + 1)] # by the names of the diseases
  
  if ( !is.null(fname) ) unlink(fname)
  enrichments.foi <- vector(mode="list", length=length(pval.foi)) # List to hold enrichment results
  # View most significant cell lines
  for(d in 1:length(pval.foi)){ # Go through each disease
    cells.foi.tmp <- pval.foi[[d]][ pval.foi[[d]] < pval]
    stats.foi.tmp <- stats.foi[[d]][ pval.foi[[d]] < pval]
    if(length(cells.foi.tmp) > 0) {
      cells.stats.foi <- as.data.frame(merge(as.matrix(cells.foi.tmp, ncol=1), t(as.data.frame(stats.foi.tmp)), by="row.names")) # Combine cell type-specific p-values with enrichment stats
      colnames(cells.stats.foi) <- c("CellType", "pval", "NumOfTests", "AvPvalCell", "AvPvalTot") # Name the columns
      cells.stats.foi <- left_join(cells.stats.foi, unique(gfAnnot[, c("cell", "cell_desc")]), by=c("CellType" = "cell")) # Join with cell annotations. We need unique to keep unique rows.
      cells.stats.foi[, c("AvPvalCell", "AvPvalTot")] <- mtx.untransform(cells.stats.foi[, c("AvPvalCell", "AvPvalTot")]) # Untransform average p-values
      # Formatting
      if (nrow(cells.stats.foi) > 1) cells.stats.foi <- cells.stats.foi[ order(cells.stats.foi$pval), ] # If more than 1 row, order by pval
      cells.stats.foi$pval <- formatC(cells.stats.foi$pval, format="e", digits=2)
      cells.stats.foi$AvPvalCell <- formatC(cells.stats.foi$AvPvalCell, format="e", digits=2)
      cells.stats.foi$AvPvalTot <- formatC(cells.stats.foi$AvPvalTot, format="e", digits=2)
      rownames(cells.stats.foi) <- cells.stats.foi$CellType
      cells.stats.foi$CellType <- NULL
      # Store results
      enrichments.foi[d] <- list(cells.stats.foi)
      if ( !is.null(fname) ) write.xlsx2(cells.stats.foi[ order(as.numeric(cells.stats.foi[, 2]), decreasing = FALSE), ], fname, sheetName=names(pval.foi)[d], row.names=FALSE, append=TRUE)
    } else {
      enrichments.foi[d] <- list("Nothing significant")
    }
    names(enrichments.foi)[d] <- names(pval.foi)[d]
  }
  return(enrichments.foi)
}

