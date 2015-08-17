Examples of visualization and analysis of the enrichment and epigenomic similarity results
========================================================

R.genomerunner is currently developed in [shiny branch](https://github.com/mdozmorov/R.genomerunner/tree/shiny)

Scripts for processing the enrichment analysis matrixes produced by [GenomeRunner Web](http://www.genomerunner.org) or [GenomeRunner](http://sourceforge.net/projects/genomerunner/). The scripts are developed in, and best used with [RStudio](http://www.rstudio.com/)

* `Analysis.Rmd and Analysis.html` - The main tutorial for the enrichment and regulatory similarity analyses and visualization in R.

* `Analysis_v2.Rmd and Analysis2.html` - The main tutorial with steps automated using `utils2.R` and `episimilarity.R` auxillary functions.

* `Enrichment_analysis.Rmd and Enrichment_analysis.html` - Demonstration of the ideas for visualization of the enrichment results.

* `episimilarity.R` - Functions for plotting, differential analyses, cell type enrichment analysis.

* `genomeRunner_file_formatting_functions.R` - Functions 'getV1...' and 'getV2...' that process LOG and detailed output of GenomeRunner, respectively, into standardized matrixes 'matrix_OR.txt' and 'matrix_PVAL.txt'. These matrixes of odds ratios and p-values are then utilized by the `utils2.R` and `episimilarity.R` functions.

* `genomeRunner_file_formatting_functions2.R` - Version 2 of the above functions that includes proper odds ratios handling.

* `utils.R` - Functions used by main scripts. Main function plotting tme majority of the results.

* `utils1.R` - Version 1 of the main function that treats the results from same cell-factor experiments from different institutions as replicates, and plots the most significant results. 

* `utils2.R` - Version 2 of the main 'showHeatmap' function, data load and other helper functions.

* `01_heatmap_corr.R` - Code snippets to perform enrichment- and regulatory similarity analyses clustering and visualization, as well as some exploratory tests.

