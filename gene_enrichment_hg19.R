# Load necessary packages
library(Biobase)
library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)
library(GOstats)
# Preparing environment for remapping Gene Symbols to Entrez IDs
x <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

## Perform GO or KEGG enrichment analysis of a list of gene symbols or EntrezIDs
##
## Parameters
##     selected - a character vector of gene IDs
##     all.universe - a all.universe of all IDs used in the study. Default - all IDs avaliable.
##     id - what type of ID is provided, "hgnc_symbol" or "entrezgene"
##     use - Which analysis to perform, "GO" or "KEGG"
##     ont - If "GO", which namespace to use, "MF", "BP", of "CC". Not used in "KEGG" analysis
##     fileName - save the ress into a file
##
## Returns a data frame of top 20, or less, significant enrichments
## Saves gene names for each of the enriched GO or KEGG at the end of the file
##
## Usage: 
##     selected <- Enrichment(gene_names_list, id="symbol", use="GO", , ont="MF")
##     selected <- Enrichment(gene_names_list, id="entrezgene", use="KEGG", "ress.txt")
##     if (nrow(selected) > 10) n <-10 else n <- nrow(selected)
##     grid.table(selected[1:n, ], gp=gpar(fontsize=7))
##
Enrichment <- function(selected, all.universe=NULL, id="symbol", use="GO", ont="BP", fileName=NULL)  {
    # Create a all.universe of all genes
    if (id == "symbol") {
      # Convert selected and all gene names to Entrez IDs, removing NAs
      selected <- unlist(xx)[selected]; selected <- selected[!is.na(selected)]
    } else if (id == "entrezgene"){
      selected <- selected
    } else {
      return("Wrong gene id type. Use 'symbol' or 'entrezgene'")
    }
  # Get universe of all genes
  if (is.null(all.universe)) {
    all.universe <- unique(c(selected, unlist(xx))); all.universe <- all.universe[!is.na(all.universe)] 
  } else {
    all.universe <- unique(c(selected, unlist(xx)[all.universe])); all.universe <- all.universe[!is.na(all.universe)] 
  }
  # Get GO-gene annotations
  geneList.annot <- select(org.Hs.eg.db,
                        keys = selected,
                        columns=c("ENTREZID", "SYMBOL", "GOALL", "PATH"),
                        keytype="ENTREZID")
  # Prepare parameters for the enrichment analysis
  if (use == "GO")
    {
    params <- new('GOHyperGParams', geneIds=selected, universeGeneIds=all.universe, ontology=ont, pvalueCutoff=0.05, conditional=F, testDirection='over', annotation="org.Hs.eg.db")
    # GO-gene summarization
    genes.annot <- split(geneList.annot[ !duplicated(geneList.annot[, 1:3]), 1:3], geneList.annot$GOALL[!duplicated(geneList.annot[, 1:3]) ]) 
    genes.annot <- lapply(genes.annot, function(x) data.frame(ID = x$GOALL[1],
                                                        SYMBOL = paste(x$SYMBOL, collapse = ","),
                                                        ENTREZID = paste(x$ENTREZID, collapse = ",")))
    genes.annot <- do.call("rbind", genes.annot) # resing data frame
    }
 else
   {
    params <- new('KEGGHyperGParams', geneIds=selected, universeGeneIds=all.universe, pvalueCutoff=0.05, testDirection='over', annotation="org.Hs.eg.db") 
    # Same for KEGG-gene summarization
    genes.annot <- split(geneList.annot[ !duplicated(geneList.annot[, c(1, 2, 6)]), c(1, 2, 6)], geneList.annot$PATH[!duplicated(geneList.annot[, c(1, 2, 6)]) ]) 
    genes.annot <- lapply(genes.annot, function(x) data.frame(ID = x$PATH[1],
                                                        SYMBOL = paste(x$SYMBOL, collapse = ","),
                                                        ENTREZID = paste(x$ENTREZID, collapse = ",")))
    genes.annot <- do.call("rbind", genes.annot) # resing data frame
   }
  # Enrichment analysis
  hgOver <- hyperGTest(params)
  res <- summary(hgOver)
  res <- cbind(res, p.adjust(res$Pvalue, method="BH")) # Append corrected for multiple testing p-value
  colnames(res)[length(colnames(res))] <- "p.adj"
  res <- res[res$p.adj < 0.1, ] # Subset the ress keeping FDR at 10%
  colnames(res)[1] <- "ID" # Set column name to merge by to ID, instead of GO- or KEGG specific
  # If genes.annot is empty, skip joining
  if ( !is.null(genes.annot)) res <- left_join(res, genes.annot, by=c("ID" = "ID")) 
  # Save the ress
  if (!is.null(fileName)) {
    write.table(res, fileName, sep="\t", row.names=F, quote=F)
  }
#   # Return only top 20 or less ress
#   ifelse(nrow(res) > 20, n <- 20, n <-nrow(res)) # Save top 20 or less enrichment ress
#   return(res[1:n, ])
  return(res)
}

# Extract and append multiple values embedded in rows, e.g. ABC11 /// BCD22 will be split into
# two separate entries, ABC11 and BCD22
#
# data: data.frame
# col: column name containing embedded values
# sep: regular expselectedsion to split column by
#
# df <- data.frame(key = c("a", "a;b", "a;b;c"), val = 1:3)
# unembed(df, "key", ";")

unembed <- function(data, col, sep, ...) {
  
  stopifnot(is.data.frame(data))
  col_i <- which(names(data) == col)
  
  data[[col]] <- as.character(data[[col]])
  pieces <- strsplit(data[[col]], sep)
  ns <- vapply(pieces, length, integer(1))
  
  #   structure(data.frame(unlist(pieces), 
  #                        data[rep(seq_along(ns), ns), -col_i]), 
  #                        names = c(col, names(data)[-col_i]))
  data.unembed <- cbind(unlist(pieces), data[rep(seq_along(ns), ns), -col_i])
  names(data.unembed) <- c(col, names(data)[-col_i])
  return(data.unembed)
}