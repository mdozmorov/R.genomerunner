Examples of visualization and analysis of the enrichment and epigenomic similarity results
========================================================

Scripts for processing the enrichment analysis matrixes produced by [GenomeRunner Web](http://www.genomerunner.org) or [GenomeRunner](http://sourceforge.net/projects/genomerunner/). The scripts are developed in, and best used with [RStudio](http://www.rstudio.com/)

* `Analysis.Rmd and Analysis.html` - The main tutorial for the enrichment and regulatory similarity analyses and visualization in R.

* `Enrichment_analysis.Rmd and Enrichment_analysis.html` - Demonstration of the ideas for visualization of the enrichment results.

* `utils.R` - Functions used by main scripts. Main function plotting tme majority of the results

* `utils1.R` - New version of the main function that treats the results from same cell-factor experiments from different institutions as replicates, and plots the most significant results. 

* `01_heatmap_corr.R` - Code snippets to perform enrichment- and regulatory similarity analyses clustering and visualization, as well as some exploratory tests.

