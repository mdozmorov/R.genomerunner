Examples of visualization and analysis of the enrichment and epigenomic similarity results
========================================================

Scripts for processing the enrichment analysis matrixes produced by [GenomeRunner Web](http://www.genomerunner.org) or [GenomeRunner](http://sourceforge.net/projects/genomerunner/). The scripts are developed in, and best used with [RStudio](http://www.rstudio.com/)

* `Analysis.Rmd and Analysis.html` - The main tutorial for the enrichment and regulatory similarity analyses and visualization in R.

* `Analysis_v2.Rmd and Analysis2.html` - The main tutorial with steps automated using `utils2.R` and `episimilarity.R` auxillary functions.

* `Enrichment_analysis.Rmd and Enrichment_analysis.html` - Demonstration of the ideas for visualization of the enrichment results.

* `episimilarity.R` - Functions for plotting, differential analyses, cell type enrichment analysis.

* `genomeRunner_file_formatting_functions.R` - Functions 'getV1...' and 'getV2...' that process LOG and detailed output of GenomeRunner, respectively, into standardized matrixes 'matrix_OR.txt' and 'matrix_PVAL.txt'. These matrixes of odds ratios and p-values are then utilized by the `utils2.R` and `episimilarity.R` functions.

* `genomeRunner_file_formatting_functions2.R` - Version 2 of the above functions that includes proper odds ratios handling.

* `gfAnnot.R` - semi-manual way of creating annotation RData object

* `stats_tests.R` - correlate GR results obtained with different statistical tests

* `utils.R` - Functions used by main scripts. Main function plotting tme majority of the results.

* `utils1.R` - Version 1 of the main function that treats the results from same cell-factor experiments from different institutions as replicates, and plots the most significant results. 

* `utils2.R` - Version 2 of the main 'showHeatmap' function, data load and other helper functions.

* `01_heatmap_corr.R` - Code snippets to perform enrichment- and regulatory similarity analyses clustering and visualization, as well as some exploratory tests.

* `jaccard` folder - BedTools script to calculate a matrix of Jaccard coefficients on all BED files in the folder

Instead of all available cell lines, subsets of tissue-specific cell lines can be used.

* Blood cell lines, any karyotype: c("K562", "Gm12878", "Gm12891", "Gm12892", "Gm06990", "Nb4", "Hl60", "Cd20", "Th1", "Gm12865", "Jurkat", "Dnd41", "Gm12864", "Th2", "Gm19239", "Cd20ro01778", "Cmk", "Gm19240", "Gm12875", "Gm12873", "Gm12872")
* Blood cell lines, normal karyotype: c("Gm12878", "Cd20ro01778", "Cd20", "Cd20ro01794")
* Blood vessel cell lines: c("Huvec", "Hpaf", "Aoaf", "Hbmec", "Hmvecdblad", "Aosmc")
* Blood and T-cell, Roadmap mnemonics: c("E033", "E034", "E037", "E038", "E039", "E040", "E041", "E042", "E043", "E044", "E045", "E047", "E048", "E062")
* HSC and B-cells, Roadmap mnemonics: c("E029", "E030", "E031", "E032", "E035", "E036", "E046", "E050", "E051")
* Brain cell lines, Roadmap mnemonics: c("E067", "E068", "E069", "E070", "E071", "E072", "E073", "E074", "E081", "E082", "E125", "E053", "E054")
* Brain cell lines, any karyotype: c("Sknsh", "Nha", "Pfsk1", "Sknmc", "Be2c", "U87", "Hah", "Gliobla", "M059j")
* Brain cell lines, normal karyotype: c("Nha", "Hah", "Bcbrainh11058n")
* Breast cell lines: c("Mcf7", "Hmec", "T47d", "Mcf10aes")
* Breast cell lines: c("E119", "E027", "E028")
* Embryonic cell lines: c("H1hesc", "Hsmm", "H1neurons", "H7es", "H9es")
* Lung cell lines: c("A549", "Imr90", "Hpae", "Hpaf", "Hpf", "Nhbe", "Nhlf", "Saec", "Wi38")

