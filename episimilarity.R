
#' Replacing 0 odds ratios
#' 
#' A function to replate odds ratios equal to zero by a minimum, or a user-specified value
#' 
#' @param mtx a matrix of odds ratios
#' @param replaceby a value to replace zero odds ratios. Default - minimum from the matrix
#'
#' @return a matrix with zeros replaced
#' @export
#' @examples
#' mtx.tumor.cor <- mtx.tumor.cor %>% set.min(replaceby=0.01)
##
set.min <- function(mtx, replaceby=min(mtx[mtx != 0])) {
  mtx[mtx == 0] <- replaceby
  return(mtx)
}

#' Heatmap plotting wrapper
#' 
#' A function to plot clustered enrichment heatmap
#' 
#' @param mtx a square correlation matrix
#'
#' @return plot clustered heatmap
#' @return complete heatmap.2 object
#' @export
#' @examples
#' mtx.tumor.cor <- mtx.tumor.cor %>% mtx.plot
##
mtx.plot <- function(mtx, dist.method="euclidean", hclust.method="ward.D2", SideColors) {
  par(oma=c(5,0,0,5), mar=c(10, 4.1, 4.1, 5), cex.main=0.65) # Adjust margins
  my.breaks <- seq(min(mtx[mtx!=min(mtx)]), max(mtx[mtx!=max(mtx)]), length.out=(2*granularity + 1))
  h <- heatmap.2(as.matrix(mtx), trace="none", density.info="none", col=color, distfun=function(x){dist(x, method=dist.method)}, hclustfun=function(x){hclust(x, method=hclust.method)}, cexRow=1, cexCol=1, breaks=my.breaks, main="Regulatory similarity clustering", RowSideColors=SideColors, ColSideColors=SideColors)  
}

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

#' Defines epigenetic marks differentially enriched in each group
#' 
#' @param mtx an episimilarity matrix 
#' @param clust clustering definition from 'mtx.clusters'
#' @param label a label to be appended to a file name. Default - "tumor"
#' @param cutoff.pval p-value cutoff to use when testing for differential enrichment. Default - 0.1
#' @param cutoff.adjust a method to correct p-values for multiple testing
#'
#' @return prints a summary of the counts of differentially enriched marks
#' @return saves the differentially enriched marks in a file
#' @export
#' @examples
#' mtx.degfs(mtx.tumor[, tumor.clust$eset.labels] %>% mtx.transform.p2z %>% normalizeQuantiles , tumor.clust, label="tumor_gfs")
##
mtx.degfs <- function(mtx, clust, label=NULL, cutoff.pval=0.1, cutoff.adjust="fdr") {
  # Limma on clusters
  eset<-new("ExpressionSet", exprs=(as.matrix(mtx)))
  # Make model matrix
  design<-model.matrix(~ 0+factor(clust$eset.groups))
  colnames(design)<-paste("c", unique(clust$eset.groups), sep="")
  # Create an empty square matrix to hold counts of DEGs
  degs.matrix<-matrix(0, length(unique(clust$eset.groups)), length(unique(clust$eset.groups)))
  colnames(degs.matrix)<-paste("c", unique(clust$eset.groups), sep="")
  rownames(degs.matrix)<-paste("c", unique(clust$eset.groups), sep="") 
  unlink(paste("results/degfs", label, ".xlsx", sep="_"))
  for(i in 1:length(colnames(design))){ 
    for(j in 1:length(colnames(design))){
      # Test only unique pairs of clusters
      if (i < j) {
        degs <- apply(exprs(eset), 1, function(x) t.test(x[design[, i] == 1], x[design[, j] == 1])$p.value)
        degs <- p.adjust(degs, method=cutoff.adjust)
        degs <- degs[degs < cutoff.pval]
        if(sum(degs < cutoff.pval) > 0) {
          # Average values in clusters i and j
          i.av<-2^rowMeans(exprs(eset)[names(degs), design[, i] == 1, drop=FALSE]) # Anti-log2 transform mean odds ratios
          j.av<-2^rowMeans(exprs(eset)[names(degs), design[, j] == 1, drop=FALSE])
          
          # Merge and convert the values
          degs.pvals <- as.matrix(cbind(degs, i.av, j.av)) 
          colnames(degs.pvals) <- c("adj.p.val", colnames(design)[i], colnames(design)[j])
          degs.pvals <- degs.pvals[order(degs.pvals[, "adj.p.val"]), , drop=FALSE]
          print(paste(colnames(design)[i], "vs.", colnames(design)[j], ", number of degs significant at adj.p.val <", cutoff.pval, ":", nrow(degs.pvals)))
          
          # Keep the number of DEGs in the matrix
          degs.matrix[i, j] <- nrow(degs.pvals)
          degs.table <- merge(degs.pvals, gfAnnot, by.x="row.names", by.y="name", all.x=TRUE, sort=FALSE) # Merge with the descriptions
          if (!is.null(label)) {
            # Format columns
            degs.table[, 2] <- formatC(degs.table[, 2], format="e", digits=2)
            degs.table[, 3] <- formatC(degs.table[, 3], format="f", digits=3)
            degs.table[, 4] <- formatC(degs.table[, 4], format="f", digits=3)
            write.xlsx2(degs.table, paste("results/degfs", label, ".xlsx", sep="_"), sheetName=paste(colnames(design)[i], "vs", colnames(design)[j], sep="_"), row.names=FALSE, append=TRUE)
          }
        }
      } 
    }
  }
  print("Counts of differential regulatory elements")
  pander(degs.matrix)
  return(degs.matrix)
}

#' Defines genes differentially expressed in each group
#' 
#' @param mtx a matrix of gene expression
#' @param clust clustering definition from 'mtx.clusters'
#' @param label a label to be appended to a file name. Default - "tumor"
#' @param cutoff.pval p-value cutoff to use when testing for differential enrichment. Default - 0.1
#' @param cutoff.adjust a method to correct p-values for multiple testing
#'
#' @return prints a summary of the counts of differentially expressed genes
#' @return saves the differentially expressed genes in a file
#' @export
#' @examples
#' mtx.degs(as.matrix(rnaseq[, normal.clust$eset.labels] + 1) %>% normalizeQuantiles %>% log2 %>% varFilter, normal.clust, label="normal_degs", cutoff.pval=0.1, cutoff.adjust="fdr")
##
mtx.degs <- function(mtx, clust, label=NULL, cutoff.pval=0.1, cutoff.adjust="fdr") {
  # Limma on clusters
  eset<-new("ExpressionSet", exprs=(as.matrix(mtx)))
  # Make model matrix
  design<-model.matrix(~ 0+factor(clust$eset.groups)) 
  colnames(design)<-paste("c", unique(clust$eset.groups), sep="")
  # Create an empty square matrix to hold counts of DEGs
  degs.matrix<-matrix(0, length(unique(clust$eset.groups)), length(unique(clust$eset.groups)))
  colnames(degs.matrix)<-paste("c", unique(clust$eset.groups), sep="")
  rownames(degs.matrix)<-paste("c", unique(clust$eset.groups), sep="") 
  unlink(paste("results/degs", label, ".txt", sep="_"))
  for(i in 1:length(colnames(design))){ 
    for(j in 1:length(colnames(design))){
      # Test only unique pairs of clusters
      if (i < j) {
        # Contrasts between two clusters
        contrast.matrix<-makeContrasts(contrasts=paste(colnames(design)[i], colnames(design)[j], sep="-"), levels=design)
        fit <- lmFit(eset, design) 
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        degs<-topTable(fit2, number=nrow(exprs(eset)), adjust.method=cutoff.adjust, p.value=cutoff.pval)
        if(nrow(degs) > 0) {
          degs <- degs[order(degs$adj.P.Val, decreasing=F), ] # Order them by the significance of the differences
          # Keep the number of DEGs in the matrix
          degs.matrix[i, j] <- nrow(degs)
          if(!is.null(label)) {
            write.table(degs, paste("results/degs", label, ".txt", sep="_"), sep="\t", quote=F,  col.names=NA, append=TRUE)
          }
        }
      } 
    }
  }
  print("Counts of differentially expressed genes")
  pander(degs.matrix)
  return(degs.matrix)
}

#' Differentia clinical parameters
#' 
#' A function to detect differential proportions of clinical data
#' 
#' @param mtx a matrix of clinical annotations, rows are patient names, columns are categorical clinical parameters
#' @param clust a list of 'eset.labels' and 'eset.groups' defining the clusters, obtained from 'mtx.clusters'
#' $param label a label to be appended to a file name. Default - "tumor"
#'
#' @return print number of differential clinical parameters in each comparison
#' @return save the differential clinical parameters into one file
#' @export
#' @examples
#' mtx.clinparam(mtx.clin, normal.clust, label="normal")
##
mtx.clinparam <- function(mtx, clust, label=NULL) {
  eset<-new("ExpressionSet", exprs=as.matrix(t(mtx[row.names(mtx) %in% clust$eset.labels, ])))
  # Make model matrix
  design<-model.matrix(~ 0+factor(clust$eset.groups)) 
  colnames(design)<-paste("c", unique(clust$eset.groups), sep="")
  # Create an empty square matrix to hold counts of DEGs
  degs.matrix<-matrix(0, length(unique(clust$eset.groups)), length(unique(clust$eset.groups)))
  colnames(degs.matrix)<-paste("c", seq(1,length(unique(clust$eset.groups))), sep="")
  rownames(degs.matrix)<-paste("c", seq(1, length(unique(clust$eset.groups))), sep="") 
  # Tweak p-value and log2 fold change cutoffs, and adjustment for multiple testing
  cutoff.pval <- 0.05
  degs <- list() # List that holds p-values for each cluster vs. cluster comparison
  ctns <- list() # List that holds counts for each cluster vs. cluster comparison
  for(i in 1:length(colnames(design))){ 
    for(j in 1:length(colnames(design))){
      # Test only unique pairs of clusters
      if (i < j) {
        d <- vector(mode="character") # List that holds p-values for current cluster vs. cluster comparison
        ct <- vector(mode="character") # List that holds counts for current cluster vs. cluster comparison
        # Test each clinical parameter
        for (z in 1:nrow(exprs(eset))) { 
          refs <- names(sort(table(exprs(eset)[z, (design[, i] == 1 | design[, j] == 1)]), decreasing=T)) # Levels of the clin. par.
          drefs <- vector(mode="numeric", length=length(refs)) # Vector to hold p-values for each level
          ctrefs <- vector(mode="character", length=length(refs)) # Vector to hold counts for each level
          # Proportion test for each level in the current clinical parameter
          for (r in 1:length(refs)) { 
            count1 <- sum(exprs(eset)[z, design[, i] == 1] == refs[r], na.rm=T)
            total1 <- sum(!is.na(exprs(eset)[z, design[, i] == 1]))
            count2 <- sum(exprs(eset)[z, design[, j] == 1] == refs[r], na.rm=T)
            total2 <- sum(!is.na(exprs(eset)[z, design[, j] == 1]))
            # Proportion test. http://ww2.coastal.edu/kingw/statistics/R-tutorials/proport.html
            drefs[r] <- prop.test(x=c(count1, count2), n=c(total1, total2), alternative="two.sided", conf.level=.99)$p.value
            ctrefs[r] <- paste(count1, "/", total1, " vs. ", count2, "/", total2, sep="") # Store counts for i cluster only
          }
          names(drefs) <- paste(rownames(exprs(eset))[z], refs, sep=".") # Which clin. par. was tested, and which level
          names(ctrefs) <- paste(rownames(exprs(eset))[z], refs, sep=".")
          d <- c(d, drefs) # Append results of current level analysis to the results of current clin. par. analysis
          ct <- c(ct, ctrefs)
        }
        # Keep the p-values
        degs <- c(degs, list(d)) # Append to list of p-values
        names(degs)[length(degs)] <- paste(colnames(design)[i], "vs", colnames(design)[j], sep="_") # Name the list
        # Keep the number of counts
        ctns <- c(ctns, list(ct)) # Append to list of counts
        names(ctns)[length(ctns)] <- paste(colnames(design)[i], "vs", colnames(design)[j], sep="_") # Name the list
      }
    }
  }
  # Save the results
  if (!is.null(label)) {
  unlink(paste("results/clinparams", label, ".txt", sep="_"))
    for (comparison in 1:length(degs)) {
      if( length(degs[[comparison]][ degs[[comparison]] < cutoff.pval]) > 0 ) {
        tmp <- cbind(degs[[comparison]][ degs[[comparison]] < cutoff.pval ],
                     ctns[[comparison]][ degs[[comparison]] < cutoff.pval ])
        colnames(tmp) <- c(names(degs)[ comparison ], names(ctns)[ comparison ])
        write.table(tmp, paste("results/clinparams", label, ".txt", sep="_"), sep="\t", quote=F, col.names=NA, append=TRUE)    
      }
    }
  }
  return(sapply(degs, function(x) length(x[ x < cutoff.pval ]))) # How many clin. par. are significant in each comparison
}

#' Randomize a matrix
#' 
#' A function to randomize a matrix using different methods
#' 
#' @param mtx a matrix of numerical values
#' @param randomize a method to randomize the matrix. 
#'        "row" (default) - replace each row with numbers sampled from a normal distribution with mean and SD of the original row. Make sure the original distribution is normal.
#'        "col" - replace each row with numbers sampled from a normal distribution with mean and SD of the original row. Make sure the original distribution is normal.
#'        "mix" - reshuffles all the numbers in the original matrix.
#'        "rnd" - replaces the entire matrix with numbers sampled from a normal distribution with mean and SD of the whole matrix.
#' @return a matrix of the same dimensions, but with randomized numbers
#' @export
#' @examples
#' mtx.rand(mtx.degs, randomize="mix")
##
mtx.rand <- function(mtx, randomize="row") {
  mtx.rnd <- matrix(NA, nrow=nrow(mtx), ncol=ncol(mtx)) # A matrix to hold random numbers
  colnames(mtx.rnd) <- colnames(mtx)
  rownames(mtx.rnd) <- rownames(mtx)
  row.sd <- rowSds(as.matrix(mtx)); col.sd <- rowSds(t(as.matrix(mtx)))
  row.mean <- rowMeans(as.matrix(mtx)); col.mean <- colMeans(as.matrix(mtx))
  if(randomize == "row") {
    for(i in 1:nrow(mtx.rnd)) {
      mtx.rnd[i, ] <- rnorm(ncol(mtx.rnd), row.mean[i], row.sd[i])
    }
  } else if(randomize == "col") {
    for (i in 1:ncol(mtx.rnd)) {
      mtx.rnd[, i] <- rnorm(nrow(mtx.rnd), col.mean[i], col.sd[i])
    }
  } else if(randomize == "mix") {
    mtx.rnd <- melt(as.matrix(mtx))
    mtx.rnd$value <- sample(mtx.rnd$value)
    class(mtx.rnd$value) <- "numeric"
    mtx.rnd <- dcast(mtx.rnd, Var1 ~ Var2, value.var="value", mean)
    rownames(mtx.rnd) <- mtx.rnd$Var1
    mtx.rnd$Var1 <- NULL
  } else if(randomize == "rnd") {
    mtx.rnd <- melt(as.matrix(mtx))
    class(mtx.rnd$value) <- "numeric"
    mtx.rnd$value <- rnorm(nrow(mtx.rnd), mean(mtx.rnd$value), sd(mtx.rnd$value))
    mtx.rnd <- dcast(mtx.rnd, Var1 ~ Var2, value.var="value", mean)
    rownames(mtx.rnd) <- mtx.rnd$Var1
    mtx.rnd$Var1 <- NULL
  }
  return(mtx.rnd)
}

#' Cell type-specific enrichment analysis
#' 
#' A function to perform cell type-specific enrichment analysis on a matrix of the enrichment results
#' 
#' @param mtx a matrix of the enrichment results. The cell type-specific enrichment analysis is performed on each column.
#' @param fname a path to a filename to save the results. Should be with "xlsx" extension. The results for each column are saved in separate worksheets.
#' @return Nothing, saves the results into the file
#' @export
#' @examples
#' mtx.cellspecific(mtx, fname="results/cellspecific_Roadmap.xlsx")
##
mtx.cellspecific <- function(mtx, fname) {
  n.diseases <- ncol(mtx) # Total number of diseases to calculate the cell type-specific p-values
  # Prepare the matrix for merging with GF annotations
  mtx <- as.data.frame(cbind(GF=rownames(mtx), mtx))
  mtx <- left_join(mtx, gfAnnot[, c("name", "cell", "factor", "description")], by=c("GF" = "name")) 
  
  pval <- 0.01 # P-value cutoff above which count the enrichments significant
  cells <- unique(mtx$cell) # All cell types
  # Global counts
  tot.tests <- nrow(mtx) # Total number of enrichment analyses
  cells.tests <- vector(mode="numeric", length=length(cells)) # Number of analyses per cell type
  for(c in 1:length(cells)) {
    cells.tests[c] <- length(mtx$cell[ mtx$cell == cells[c]]) # Number of analyses per cell type
  }
  
  # Disease-specific counts
  pval.disease <- list() # for p-values
  c2x2.disease <- list() # for 2x2 tables
  for(i in 2:(n.diseases + 1)){ # Disease positions are now shifted by 1
    tot.sig <- length(mtx[ 1/10^abs(as.numeric(mtx[, i])) < pval, i]) # Disease-specific total number of significant results
    cells.sig <- vector(mode="numeric", length=length(cells)) # Disease-specific & cell type-specific number of significant results
    for(c in 1:length(cells)) {
      mtx.sel <- as.numeric(mtx[ mtx$cell == cells[c], i ]) # Disease- and cell type-specific vector
      cells.sig[c] <- length( mtx.sel[ 1/10^abs(mtx.sel) < pval ]  ) # How many are significant
    }
    pval.disease.cell <- vector(mode="numeric", length=length(cells)) # A vector to store disease- and cell type-specific p-values
    c2x2.disease.cell <- list() # A list to store disease- and cell type-specific 2x2 tables
    for(c in 1:length(cells)) {
      cells.c2x2 <- matrix(c(cells.sig[c], cells.tests[c] - cells.sig[c], 
                             tot.sig - cells.sig[c], tot.tests - tot.sig - (cells.tests[c] - cells.sig[c])),
                           dimnames=list(c("cell", "not cell"), c("sig", "not sig")),
                           ncol=2) # 2x2 contingency table
      cells.test <- fisher.test(cells.c2x2)
      if(cells.test$estimate > 1){
        pval.disease.cell[c] <- cells.test$p.value # Store the enrichment p-values
      } else {
        pval.disease.cell[c] <- 1 # If depletion, it's not interesting
      }
      c2x2.disease.cell <- c( c2x2.disease.cell, list(c(cells.test$estimate, cells.sig[c], cells.tests[c], tot.sig, tot.tests))) # Store the numbers to construct 2x2 tables
    }
    names(pval.disease.cell) <- cells # Name the collected vectors
    names(c2x2.disease.cell) <- cells # as cell names
    pval.disease <- c(pval.disease, list(pval.disease.cell)) # Store them
    c2x2.disease <- c(c2x2.disease, list(c2x2.disease.cell)) # in disease-specific lists
  }
  names(pval.disease) <- colnames(mtx)[2:(n.diseases + 1)] # Name the disease-specific lists
  names(c2x2.disease) <- colnames(mtx)[2:(n.diseases + 1)] # by the names of the diseases
  
  unlink(fname)
  # View most significant cell lines
  for(d in 1:length(pval.disease)){ # Go through each disease
    print(names(pval.disease)[d])
    print(d)
    cells.disease <- pval.disease[[d]][ pval.disease[[d]] < pval]
    stats.disease <- c2x2.disease[[d]][ pval.disease[[d]] < pval]
    if(length(cells.disease) > 0) {
      # cells.stats.disease <- cbind(cells.disease, ldply(stats.disease, rbind))
      cells.stats.disease <- as.data.frame(merge(as.matrix(cells.disease, ncol=1), t(as.data.frame(stats.disease)), by="row.names"))
      rownames(cells.stats.disease) <- cells.stats.disease$Row.names
      cells.stats.disease$Row.names <- NULL
      colnames(cells.stats.disease) <- c(names(pval.disease)[d], "OR", "cell.sig", "cell.tot", "all.sig", "all.tot")
      cells.stats.disease <- merge(cells.stats.disease, cellAnnot, by.x="row.names", by.y="cell")
      class(cells.stats.disease$description) <- "character"
      cells.stats.disease[, 2] <- formatC(cells.stats.disease[, 2], format="e", digits=2)
      cells.stats.disease[, 3] <- formatC(cells.stats.disease[, 3], format="f", digits=2)
      write.xlsx2(cells.stats.disease[ order(as.numeric(cells.stats.disease[, 3]), decreasing = TRUE), ], fname, sheetName=names(pval.disease)[d], row.names=FALSE, append=TRUE)
    } else {
      #write.xlsx2(names(pval.disease)[d], fname,  sheetName=names(pval.disease)[d], row.names=FALSE, append=TRUE)
      write.xlsx2("No cell type-specific enrichments", fname, sheetName=names(pval.disease)[d], row.names=FALSE, append=TRUE)
    }
  }
}