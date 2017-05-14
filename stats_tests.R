# Matrixes of p-values and odds ratios
# Chi-square results
mtx.c.p <- read.table("data/test_chisquare/matrix_PVAL.txt", sep = "\t", stringsAsFactors = F)
mtx.c.o <- read.table("data/test_chisquare/matrix_OR.txt", sep = "\t", stringsAsFactors = F)
# Binomial results
mtx.b.p <- read.table("data/test_binomial/matrix_PVAL.txt", sep = "\t", stringsAsFactors = F)
mtx.b.o <- read.table("data/test_binomial/matrix_OR.txt", sep = "\t", stringsAsFactors = F)
# Monte Carlo results
mtx.m.p <- read.table("data/test_montecarlo/matrix_PVAL.txt", sep = "\t", stringsAsFactors = F)
mtx.m.o <- read.table("data/test_montecarlo/matrix_OR.txt", sep = "\t", stringsAsFactors = F)

# Testing correlations
cor(mtx.c.p, mtx.b.p) %>% diag %>% as.vector
cor(mtx.c.p[1:nrow(mtx.m.p), ], mtx.m.p) %>% diag %>% as.vector
