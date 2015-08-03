#' Define clusters
#' 
#' A function to define clusters by cutting clustering dendrogram produced by 'mtx.plot'. 
#' 
#' @param mtx.colDendrogram a colDendrogram part of heatmap.2 object
#' @param height a threshold to cut the dendrogram. Estimate from calling 'mtx.tumor.cor$colDendrogram %>% plot(horiz=T)'
#' @param minmembers a minimum number of samples to be considered as a cluster
#' @param label a label to be appended to a file name. Default - "tumor"
#'
#' @return a list with 'eset.labels' and 'eset.groups' slots defining clustered labels/groups from the dendrogram cutting
#' @return save the clustering order in a file
#' @export
#' @examples
#' tumor.clust <- mtx.tumor.cor$colDendrogram %>% mtx.clusters(height=1000, label="tumor")
##
mtx.clusters <- function(mtx.colDendrogram, height=10, minmembers=3, label=NULL) {
  c <- cut(mtx.colDendrogram, h=height) 
  # Check the number of clusters, and the number of members. Output the results into a file
  unlink(paste("results/clust_", label, ".txt", sep=""))
  for (i in 1:length(c$lower)) {
    cat(paste("Cluster", formatC(i, width=2, flag="0"), sep=""), "has", formatC(attr(c$lower[[i]], "members"), width=3), "members", "\n")
    if (!is.null(label)) {
      write.table(paste(i, t(labels(c$lower[[i]])), sep="\t"), paste("results/clust_", label, ".txt", sep=""), sep="\t", quote=F,  col.names=F, row.names=F, append=T)
    }
  }
  # Define Groups
  eset.labels <- character() # Empty vector to hold cluster labels
  eset.groups <- numeric() # Empty vector to hold cluster groups
  for (i in 1:length(c$lower)) { # Go through each cluster
    if (attr(c$lower[[i]], "members") > minmembers) { # If the number of members is more than a minimum
      eset.labels<-append(eset.labels, labels(c$lower[[i]]))
      eset.groups<-append(eset.groups, rep(i, length(labels(c$lower[[i]]))))
    }
  }
  return(list(eset.labels=eset.labels, eset.groups=eset.groups))
}