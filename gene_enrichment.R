# Load necessary packages
library(Biobase)
library(reactome.db)
library(KEGG.db)
#library(arrayQualityMetrics)
#library(limma)
#library(genefilter)
#library(sva)
#library(WGCNA)
#library(gplots)
#library(RColorBrewer)
#library(tspair)
library(gridExtra)
library(GO.db)
library(GOstats)
library(Hmisc)

# Human annotations
library(org.Hs.eg.db)
# Preparing environment for remapping Gene Symbols to Entrez IDs
x <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

## Perform GO or KEGG enrichment analysis of a list of gene symbols or EntrezIDs
##
## Parameters
##     res - a character vector of gene IDs
##     id - what type of ID is provided, "hgnc_symbol" or "entrezgene"
##     use - Which analysis to perform, "GO" or "KEGG"
##     ont - If "GO", which namespace to use, "MF", "BP", of "CC". Not used in "KEGG" analysis
##     fileName - save the results into a file
##
## Returns a data frame of top 20, or less, significant enrichments
## Saves gene names for each of the enriched GO or KEGG at the end of the file
##
## Usage: 
##     res <- Enrichment(gene_names_list, id="symbol", use="GO", , ont="MF")
##     res <- Enrichment(gene_names_list, id="entrezgene", use="KEGG", "results.txt")
##     if (nrow(res) > 10) n <-10 else n <- nrow(res)
##     grid.table(res[1:n, ], gp=gpar(fontsize=7))
##
Enrichment <- function(res, id="symbol", use="GO", ont="BP", fileName=NULL)
{
    if (id == "symbol"){
    # Convert selected and all gene names to Entrez IDs, removing NAs
    sel.entrez <- unlist(xx)[res]; sel.entrez <- sel.entrez[!is.na(sel.entrez)]
    all.entrez <- unlist(xx); all.entrez <- all.entrez[!is.na(all.entrez)]    
  } else if (id == "entrezgene"){
    sel.entrez <- res
    all.entrez <- unique(c(sel.entrez, unlist(xx)))
  } else {
    return("Wrong gene id type. Use 'hgnc_symbol' or 'entrezgene'")
  }
  # Prepare parameters for the enrichment analysis
  if (use == "GO")
    {
    params <- new('GOHyperGParams', geneIds=sel.entrez, universeGeneIds=all.entrez, ontology=ont,
 pvalueCutoff=0.05, conditional=F, testDirection='over', annotation="org.Hs.eg.db")
    }
 else
   {
    params <- new('KEGGHyperGParams', geneIds=sel.entrez, universeGeneIds=all.entrez, pvalueCutoff=0.05, testDirection='over', annotation="org.Hs.eg.db") 
   }
  hgOver <- hyperGTest(params)
  result <- summary(hgOver)
  result <- cbind(result, p.adjust(result$Pvalue, method="BH")) # Append corrected for multiple testing p-value
  colnames(result)[length(colnames(result))] <- "p.adj"
  result <- result[result$p.adj < 0.1, ] # Subset the results keeping FDR at 10%
  if (!is.null(fileName)) {
    write.table(result, fileName, sep="\t", row.names=F)
    # In addition to saving the enrichment results
    # save genes in each GO category
    if (id == "entrezgene") { keytype <- "ENTREZID" }
    if (id == "symbol") { keytype <- "SYMBOL" }
    if (use == "GO" & nrow(result) > 0) {
        geneList.go <- select(org.Hs.eg.db,
                           keys = res,
                           columns=c("ENTREZID","SYMBOL","GOALL"),
                           keytype=keytype)
        for (i in 1:nrow(result)) {
          write.table(paste(result[i, 1], paste(sort(unique(geneList.go$SYMBOL[ geneList.go$GOALL == result[i, 1]])), collapse=","), collapse="\t"), fileName, sep="\t", col.names=F, row.names=F, append=T)
        }
    }
    if (use == "KEGG" & nrow(result) > 0) {
        geneList.go <- select(org.Hs.eg.db,
                           keys = res,
                           columns=c("ENTREZID","SYMBOL","PATH"),
                           keytype=keytype)
        for (i in 1:nrow(result)) {
          write.table(paste(result[i, 1], paste(sort(unique(geneList.go$SYMBOL[ geneList.go$PATH == result[i, 1]])), collapse=","), collapse="\t"), fileName, sep="\t", col.names=F, row.names=F, append=T)
        }
    }
  }
  ifelse(nrow(result) > 20, n <- 20, n <-nrow(result)) # Save top 20 or less enrichment results
  return(result[1:n, ])
}