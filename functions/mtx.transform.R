## Convert a matrix of raw p-values (with "-" indicating depletion) into -log10-transformed
mtx.transform<-function(x){
  tmp<- -log10(abs(x)) # -log10 transformation without sign
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (x[i, j] < 0 & !is.na(x[i, j])) {tmp[i, j] <- -tmp[i,j]} # Add sign, if needed
    }
  }
  return(tmp)
}