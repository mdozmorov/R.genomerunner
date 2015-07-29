#' Defines epigenetic marks differentially enriched in each group
#' 
#' @param mtx the matrix of transformed p-values or odds-ratios
#' @param clust clustering definition from 'mtx.clusters'
#' @param label a label to be appended to a file name. Default - NULL, do not write into the file
#' @param cutoff.pval p-value cutoff to use when testing for differential enrichment. Default - 0.1
#' @param cutoff.adjust a method to correct p-values for multiple testing
#' @param isOR logical. If true, the transformed values are log2-untransformed (to become odds-ratios), otherwise, -log10-untransformed (to become p-values). Default - FALSE
#'
#' @return prints a summary of the counts of differentially enriched marks
#' @return saves the differentially enriched marks in a file
#' @export
#' @examples
#' mtx.degfs(mtx.tumor[, tumor.clust$eset.labels] %>% mtx.transform.p2z %>% normalizeQuantiles , tumor.clust, label="tumor_gfs")
##
mtx.degfs <- function(mtx, clust, label=NULL, cutoff.pval=0.1, cutoff.adjust="fdr", isOR=FALSE) {
  # Limma on clusters
  eset<-new("ExpressionSet", exprs=(as.matrix(mtx[ , clust$eset.labels])))
  # Make model matrix
  design<-model.matrix(~ 0+factor(clust$eset.groups))
  colnames(design)<-paste("c", unique(clust$eset.groups), sep="")
  # Create an empty square matrix to hold counts of DEGs
  degs.matrix<-matrix(0, length(unique(clust$eset.groups)), length(unique(clust$eset.groups)))
  colnames(degs.matrix)<-paste("c", unique(clust$eset.groups), sep="")
  rownames(degs.matrix)<-paste("c", unique(clust$eset.groups), sep="") 
  unlink(paste("results/degfs_", label, ".xlsx", sep=""))
  degfs.list <- list()
  for(i in 1:length(colnames(design))){ 
    for(j in 1:length(colnames(design))){
      # Test only unique pairs of clusters
      if (i < j) {
        degs <- apply(exprs(eset), 1, function(x) wilcox.test(x[design[, i] == 1], x[design[, j] == 1])$p.value)
        degs <- degs[ !is.na(degs) ] # Precaution against NA p-values, when both clusters have exactly the same numbers
        degs <- p.adjust(degs, method=cutoff.adjust)
        degs <- degs[degs < cutoff.pval]
        # Average values in clusters i and j
        if(sum(degs < cutoff.pval) > 0) {
          if( isOR == FALSE ) {
            i.av<-1/(10^rowMeans(abs(exprs(eset)[names(degs), design[, i] == 1, drop=FALSE]))) # Anti -log10 transform p-values
            j.av<-1/(10^rowMeans(abs(exprs(eset)[names(degs), design[, j] == 1, drop=FALSE])))
          } else {
            i.av<-2^rowMeans(exprs(eset)[names(degs), design[, i] == 1, drop=FALSE]) # Anti log2 transform mean odds ratios
            j.av<-2^rowMeans(exprs(eset)[names(degs), design[, j] == 1, drop=FALSE])
          }
          
          # Merge and convert the values
          degs.pvals <- as.matrix(cbind(degs, i.av, j.av)) 
          colnames(degs.pvals) <- c("adj.p.val", colnames(design)[i], colnames(design)[j])
          degs.pvals <- degs.pvals[order(degs.pvals[, "adj.p.val"]), , drop=FALSE]
          print(paste(colnames(design)[i], "vs.", colnames(design)[j], ", number of degs significant at adj.p.val <", cutoff.pval, ":", nrow(degs.pvals)))
          
          # Keep the number of DEGs in the matrix
          degs.matrix[i, j] <- nrow(degs.pvals)
          degs.table <- merge(degs.pvals, gfAnnot, by.x="row.names", by.y="file_name", all.x=TRUE, sort=FALSE) # Merge with the descriptions
          # Format columns
          degs.table[, 2] <- formatC(degs.table[, 2], format="e", digits=2)
          degs.table[, 3] <- formatC(degs.table[, 3], format="f", digits=3)
          degs.table[, 4] <- formatC(degs.table[, 4], format="f", digits=3)
          # Save the results in the list
          degfs.list <- c(degfs.list, list(degs.table))
          names(degfs.list)[length(degfs.list)] <- paste(colnames(design)[i], "vs", colnames(design)[j], sep="_") # Name the list after combination of clustering
          # Save the results in the file
          if (!is.null(label)) {
            write.xlsx2(degs.table, paste("results/degfs_", label, ".xlsx", sep=""), sheetName=paste(colnames(design)[i], "vs", colnames(design)[j], sep="_"), row.names=FALSE, append=TRUE)
          }
        }
      } 
    }
  }
  print("Counts of differential regulatory elements")
  pander(degs.matrix)
  return(degfs.list) # Return the full list of the results
}
