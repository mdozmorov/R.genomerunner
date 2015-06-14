suppressMessages(source("utils2.R")) # See the required packages there
suppressMessages(source("episimilarity.R"))
mtx <- load_gr_data("/home/lukas/gftest/matrix_OR.txt")
mtx <- scale(mtx)
mtx.cor <- rcorr(as.matrix(mtx), type="spearman")[[1]]
mtx.plot(mtx.cor)
